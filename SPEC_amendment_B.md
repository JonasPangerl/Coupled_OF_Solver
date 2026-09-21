<!-- Amendment set B to SPEC_coupledFoam.md (v2), received from the user on
     2026-09-21, stored verbatim. Binding: it overrides the sections it names.
     Implementation checklist: AMENDMENT_B_STATUS.md. -->

## Amendment set B — supersedes amendment set A entirely (A1–A7 are void)

### B0 — Scope of this amendment
Replaces Section 6.1, extends 6.3 (new 6.3.1–6.3.5), adds 7.1 trigger, adds
7.5, adds 8.3, extends 11, 13, 15, and updates the memory table in Section 1.
Everything not mentioned here is unchanged from SPEC v2.

### B1 — Section 6.1 replaced: runtime selection
```
coupled
{
    solver          blockFGMRES;     // blockFGMRES | blockGMRES | blockBiCGStab
    preconditioner  blockGAMG;       // blockDiagonal | blockGAMG
    tolerance       1e-8;            // absolute floor on ||r||/normFactor
    relTol          0.05;            // used only if adaptiveRelTol no
    adaptiveRelTol  yes;             // Eisenstat-Walker, see B2
    etaMin 1e-3; etaMax 0.5; gammaEW 0.9; alphaEW 2;
    minIter         1;
    maxIter         200;
    restart         10;              // blockFGMRES/blockGMRES; use 6 above 40M cells
    blockGAMG { ... see B3 ... }
}
```
Constraint enforced at start-up (FatalError otherwise):
`blockGAMG.cycleType K` requires `solver blockFGMRES`. Reason: the K-cycle is
a variable preconditioner; BiCGStab and standard GMRES assume a fixed one.

`blockFGMRES`: flexible GMRES (Saad, Alg. 9.6). Stores V_j and Z_j = M_j^{-1} V_j
for j = 1..restart → 2·restart block vectors. Modified Gram–Schmidt in double
accumulation (doubleReduce). Happy-breakdown guard: if h_{j+1,j} < VSMALL·‖r0‖
stop and solve the small system. Small Hessenberg least-squares via Givens
rotations in double.

### B2 — Adaptive inner tolerance (Eisenstat–Walker), Section 6.1 addition
Per outer iteration n (R_n as in 7.1), the linear solve target is
`‖r_lin‖ ≤ eta_n · ‖r_lin,0‖` with
  eta_n = gammaEW · (R_n / R_{n-1})^alphaEW
safeguards, applied in this order:
  1. if gammaEW · eta_{n-1}^alphaEW > 0.1: eta_n = max(eta_n, gammaEW · eta_{n-1}^alphaEW)
  2. eta_n = clamp(eta_n, etaMin, etaMax)
  3. during start-up (7.4, startupDone == false): eta_n = etaMax
  4. iteration n = 1: eta_1 = etaMax
The absolute `tolerance` remains a floor. `minIter 1` is mandatory (never a
zero increment). Log: `eta=<v>` in the CF| line.
Ref: Eisenstat & Walker, SIAM J. Sci. Comput. 17 (1996), Choice 2;
PTC application: Appl. Math. Comput. 266 (2015).

### B3 — Section 6.3 extended

#### 6.3.1 Agglomeration and coarsening-ratio rule
```
blockGAMG
{
    agglomerator          faceAreaPair;
    mergeLevels           2;          // MUST be >= 2 for cycleType W or K (see rule)
    nCellsInCoarsestLevel 200;
    processorAgglomerator masterCoarsest;   // native GAMGProcAgglomeration
    nProcessorsPerLevel   ... ;       // see 6.3.4
    cycleType             K;          // V | F | W | K
    kCycleThreshold       0.25;       // t in 6.3.2
    kCycleMaxSteps        2;
    smoother              blockGaussSeidel;
    nPreSweeps            1;
    nPostSweeps           2;
    nFinestSweeps         2;
    coarsestSolver        blockBiCGStab;
    coarsestTolerance     1e-3;
    coarsestMaxIter       50;
    maxOperatorComplexity 1.5;
    cacheAgglomeration    true;
    autoTune              yes;        // see 6.3.5
    tuneInterval          50;
    nPostSweepsMax        4;
}
```
**Coarsening-ratio rule (enforced, FatalError otherwise):** let
`r_l = nCells(l) / nCells(l+1)` be the measured ratio between consecutive
levels after agglomeration. For `cycleType W` or `K` every level must satisfy
`r_l >= 3.0`; if violated, the solver increases `mergeLevels` by 1 (max 4) and
re-agglomerates; if still violated, FatalError with the measured ratios.
Rationale: cycle cost per finest-level work unit is
  V: Σ r^{-l} ≈ r/(r−1);  F: Σ (l+1) r^{-l};  W: Σ 2^l r^{-l}  (diverges for r ≤ 2).
With pairwise agglomeration r = 2, so W/K cost O(L·N) unless levels are merged.

#### 6.3.2 Cycle algorithms (exact)
Notation: level 0 = finest. A_l block matrix, S_l(x,b,ν) = ν sweeps of the
smoother, R_l restriction (block sum), P_l prolongation (piecewise constant),
L = coarsest level index.
```
CYCLE(l, b) -> x                       // solves A_l x ≈ b, x starts at 0
  if l == L:
      x = coarsestSolver(A_L, b, coarsestTolerance, coarsestMaxIter); return x
  x = 0
  x = S_l(x, b, nPreSweeps)            // pre-smoothing (nFinestSweeps if l==0)
  r = b − A_l x
  b_c = R_l r
  switch cycleType:
    V: e_c = CYCLE(l+1, b_c)
    F: e_c = CYCLE_F(l+1, b_c)         // see below
    W: e_c = CYCLE(l+1, b_c); e_c += CYCLE(l+1, b_c − A_{l+1} e_c)
    K: e_c = KSTEP(l+1, b_c)           // see below
  x += P_l e_c
  x = S_l(x, b, nPostSweeps)           // post-smoothing (nFinestSweeps if l==0)
  return x

CYCLE_F(l, b): one W-style recursion at this level, then V below:
  identical to CYCLE(l,b) but with the recursive call pattern
  "e_c = CYCLE_F(l+1,b_c); e_c += CYCLE_V(l+1, b_c − A_{l+1} e_c)" where
  CYCLE_V is CYCLE forced to V. (Standard F-cycle: one extra visit per level.)

KSTEP(l, b) -> e                       // Notay–Vassilevski K-cycle, GCR(kCycleMaxSteps)
  // one or two steps of GCR on A_l e = b, preconditioned by CYCLE(l, ·) with V-recursion
  // below this level (i.e. CYCLE at level l uses cycleType K again, so the
  // adaptivity is applied at every level)
  r0 = b ; ‖r0‖ computed in double
  z1 = CYCLE(l, r0)                    // preconditioner application
  q1 = A_l z1
  a1 = <q1, r0> / <q1, q1>             // double-accumulated dots
  e  = a1 z1 ; r1 = r0 − a1 q1
  if ‖r1‖ <= kCycleThreshold · ‖r0‖ or kCycleMaxSteps == 1: return e
  z2 = CYCLE(l, r1)
  q2 = A_l z2
  // GCR orthogonalisation against q1
  q2 = q2 − (<q2,q1>/<q1,q1>) q1 ;  z2 = z2 − (<q2_orig,q1>/<q1,q1>) z1   [use the
      pre-orthogonalisation q2 in the coefficient]
  a2 = <q2, r1> / <q2, q2>
  e += a2 z2
  return e
```
Workspace per level for KSTEP: 4 block vectors (z1,q1,z2,q2) + r; allocated
once and cached. K-cycle uses **kCycleMaxSteps = 2** at every level except
the finest (the finest is driven by the outer FGMRES).
Ref: Notay & Vassilevski, *Recursive Krylov-based multigrid cycles*,
Numer. Linear Algebra Appl. 15 (2008).

#### 6.3.3 Coarsest solver
`blockBiCGStab` with blockDiagonal preconditioner, tolerance 1e-3 relative,
maxIter 50; if the coarsest level has ≤ 64 cells on the owning rank, use a
dense LU of the assembled (4·n)×(4·n) matrix instead (double, factorised once
per outer iteration when the matrix changes). Selection is automatic and
logged.

#### 6.3.4 Processor agglomeration (coarse-level latency)
Coarse levels are gathered onto fewer ranks using native
`GAMGProcAgglomeration` (`masterCoarsest` with `nProcessorsPerLevel`):
default rule — a level whose global cell count is below `procAgglomCells`
(= 20000 · nRanks... no: use 5000 cells per rank) i.e.
`nCells(l) < 5000 · nRanks` is agglomerated to `max(1, nRanks/4)` ranks, and
below `5000` cells to 1 rank. The block interface coefficients are gathered
with the same mapping as the native scalar case (`GAMGAgglomeration::
procAgglomerateLduAddressing`), with 16 coefficients per face instead of 1.
Levels and rank counts are printed once per run. This is the only sanctioned
form of "doing less on coarse levels"; skipping coarse levels is prohibited
(it removes the low-frequency correction and saves <5 % of cycle work).

#### 6.3.5 Cycle-efficiency controller (autoTune)
Observable: per outer FGMRES iteration, rho = ‖r_after_precon_step‖/‖r_before‖
is NOT directly available inside FGMRES; therefore measure rho on the
**first** preconditioner application of each linear solve as
rho = ‖b − A M^{-1} b‖ / ‖b‖ (one extra Amul per outer iteration; cost ≈ 1 %).
Every `tuneInterval` outer iterations, on the median rho of the window:
  - rho > 0.7: nPostSweeps += 1 (≤ nPostSweepsMax). If already at max and
    cycleType ∈ {V,F}: promote V→F→W. If cycleType K: nothing further
    (K-cycle adapts internally; record `kcycle-saturated` event).
  - rho < 0.3 and nPostSweeps > 1: nPostSweeps −= 1. If nPostSweeps == 1 and
    rho < 0.2 and cycleType was promoted by the controller: demote one step.
  - A change requires the condition in two consecutive windows (hysteresis).
  - Every change: log line `GAMG-tune: <what> rho_med=<v>`; JSON event.
  - Failure condition: rho_med > 0.9 for two windows at nPostSweepsMax with
    W or K → WARNING + event; on T0–T5 this is a test failure.
The controller never changes anything *during* a linear solve; changes apply
at the next outer iteration.

### B4 — Section 7.1 addition: linear-solver divergence trigger
If the linear solve exits with ‖r_lin‖ > eta_n·‖r_lin,0‖ after maxIter:
CFL ← max(CFLmin, kappa·CFL); arm nHold; repeat the outer step with unchanged
fields (counts toward maxCflCuts). After `maxLinFails` (3) consecutive
failures: abort through the 9.3 diagnostic path. Recovery afterwards is
governed by mRDM (bounded by betaMax), not by a fixed ramp.

### B5 — New Section 7.5: Anderson acceleration (optional, default off)
```
anderson { enabled no; m 4; beta 1.0; }
```
Applied to the accepted outer update after the line search:
  x_k = current (U,p) after update; g_k = x_k; w_k = ω·dx (accepted increment)
  keep last m pairs (w_j, Δw_j). Solve min_α ‖w_k − Σ_j α_j Δw_j‖ by QR of
  the m×m Gram system in double (dimension m, cost negligible), then
  x_{k+1} = x_k + β·(w_k − Σ α_j Δw_j)   (Walker & Ni form, Type II)
Safeguards (all mandatory): the extrapolated state passes the 9.3 sentinel,
else the un-extrapolated x_k is kept and the history is flushed; history is
flushed on every CFL cut, every remediation-set change, and every rollback;
the Anderson update is skipped (not flushed) if |α|_max > 10.
Memory: m pairs of (U,p) → 2·m·4 scalar fields (m = 4: 32 scalars/cell in
double = 256 B/cell ≈ 11.5 GB at 45 M cells) — therefore default off and
documented in the memory table. Benchmark: T1 and T3-SST additionally with
`enabled yes`; report the iteration and wall-clock delta (15.5).
Refs: Pollock, Rebholz et al., SIAM J. Numer. Anal. (2019); Walker & Ni,
SIAM J. Numer. Anal. 49 (2011).

### B6 — New Section 8.3: zonal factors
```
remediation.zonal
{
    enabled no;
    zones ( { cellZone wheelHousing; cflFactor 0.5; beta 0.5; } );
    patchDistance ( { patches (body); nLayers 3; cflFactor 0.7; beta 1.0; } );
}
```
Multipliers on local Δt (cflFactor) and convection blending (beta) for cells
in named cellZones or within nLayers face-neighbour layers of named patches
(layers computed once by BFS from patch faces). Multiple factors multiply.
Design note for DECISIONS.md: local pseudo-time is implicit local
under-relaxation; explicit per-patch ω factors are not implemented (they
break increment consistency of the line search). Processor boundaries need
no relaxation (block interfaces are fully implicit).

### B7 — Section 1 memory table replaced (45 M cells, worst case)
| Item | ~GB |
|---|---|
| Mesh, fields, addressing (double) | 40–45 |
| Block matrix 4×4 (single) | ~20 |
| Block-GAMG hierarchy, C_op ≤ 1.5 (single) | 20–30 |
| FGMRES restart 6: 12 block vectors (single) | ~8.6 |
| K-cycle workspace all levels (≤ 1.5× finest: 5 block vectors) | ~5 |
| Turbulence, rollback copies, reserve | ~10 |
| **Total** | **~105–120** |
| Anderson m=4 (optional) | +11.5 → not enabled above 35 M cells |
Rule: above 40 M cells set `restart 6`; below, `restart 10`.

### B8 — Section 11 keyword additions
| Keyword | Type | Default | Section |
|---|---|---|---|
| `solver` | word | blockFGMRES | B1 |
| `restart` | label | 10 | B1 |
| `adaptiveRelTol / etaMin / etaMax / gammaEW / alphaEW` | bool/scalar | yes / 1e-3 / 0.5 / 0.9 / 2 | B2 |
| `blockGAMG.cycleType / kCycleThreshold / kCycleMaxSteps` | word/scalar/label | K / 0.25 / 2 | 6.3.2 |
| `blockGAMG.mergeLevels / nCellsInCoarsestLevel` | label | 2 / 200 | 6.3.1 |
| `blockGAMG.processorAgglomerator` | word | masterCoarsest | 6.3.4 |
| `blockGAMG.coarsestMaxIter` | label | 50 | 6.3.3 |
| `blockGAMG.autoTune / tuneInterval / nPostSweepsMax` | bool/label | yes / 50 / 4 | 6.3.5 |
| `ptc.maxLinFails` | label | 3 | B4 |
| `anderson.enabled / m / beta` | bool/label/scalar | no / 4 / 1.0 | B5 |
| `remediation.zonal.enabled` | bool | no | B6 |

### B9 — Section 6.4 unit-test additions
- `Test-blockGAMG` extended: solve the block-Poisson system with cycleType
  V, F, W, K on the 128×128 cavity mesh and on motorBike-tutorial mesh
  (350 k); record iterations and wall-clock; pass: K ≤ W ≤ V in iterations;
  all converge to 1e-8; identical solution 1 vs 4 ranks to 1e-5; measured
  r_l ≥ 3 with mergeLevels 2.
- `Test-blockFGMRES`: with a deliberately variable preconditioner (random
  sweep count per application) FGMRES must converge; standard GMRES on the
  same setup is allowed to stall (documents why flexibility is required).
- `Test-procAgglom`: 16 ranks, motorBike-tutorial mesh; coarse-level rank
  counts printed match rule 6.3.4; solution identical to no-agglomeration
  run to 1e-5.

### B10 — Section 13/14/15 additions
- Benchmark matrix gains configuration (E): `cycleType V` vs default K
  (isolates the cycle benefit) on T2, T4b, T5.
- T1 and T3-SST additionally with `adaptiveRelTol no`: adaptive must not be
  > 5 % slower and must reach identical Cd/Δp within 1e-4.
- T1, T3-SST additionally with `anderson.enabled yes`.
- Report 15.7 gains: eta history, rho history with tune events, per-level
  cell/rank table, cycle-type comparison bars, Anderson on/off comparison.

---

Author's notes accompanying the amendment (verbatim summary, translated):
three points corrected against an earlier draft — the controller may no
longer switch V→W without an affordable coarsening ratio (hard rule
r_l >= 3); rho is measured explicitly on the first preconditioner
application of each solve, with its cost stated; Anderson memory is 32
scalars per cell at m = 4, hence the cell-count limit in the memory table.
Design rationale: additive-correction/aggregation AMG has a weak coarse-grid
correction, so a pure V-cycle does not scale optimally; the K-cycle
(Notay & Vassilevski 2008, basis of AGMG) switches between V- and W-like
behaviour per level and per application, based on measurement, and requires
a flexible outer Krylov method (FGMRES). Coarse levels are cheap (< 5 % of
cycle work) and must not be skipped; on many ranks their cost is MPI
latency, addressed by processor agglomeration.
