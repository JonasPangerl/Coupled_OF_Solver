# DECISIONS

Every deviation from `SPEC_coupledFoam.md`, every spec defect found, and every
interpretation of an ambiguous passage. Format: ID, status, what, why, impact.

---

## D-001 - Precision: `solveScalar` cannot be float while `scalar` is double

**Status:** open - recommendation below, user confirmation requested (2026-09-21).

**Finding.** Spec 3 requires `WM_PRECISION_OPTION=SPDP` and a Phase 0 check
`sizeof(Foam::solveScalar)==4 && sizeof(Foam::scalar)==8`. In OpenFOAM v2606
(`src/OpenFOAM/primitives/Scalar/scalar/scalarFwd.H`):

    #if defined(WM_SP)      scalar=float   solveScalar=float
    #elif defined(WM_SPDP)  scalar=float   solveScalar=double
    #elif defined(WM_DP)    scalar=double  solveScalar=double

No build option satisfies the literal criterion. SPDP is the inverse of the
spec's intent (fields in float, linear solve in double). Measured on the
installed system build: `linux64GccDPInt32Opt`, `sizeof(scalar)=8`,
`sizeof(solveScalar)=8`.

**Decision (recommended).** Keep OpenFOAM in DP. The library defines its own
storage types in `src/blockMatrix/blockScalar.H`:
`blockScalar = floatScalar` (block coefficients, Krylov/AMG vectors) and
`reduceScalar = doubleScalar` (all accumulations, MPI reductions). Where the
spec says `solveScalar` for block data, the implementation uses
`blockScalar`. Native `solveScalar` is never used for block data.
Phase 0 gate evaluates `sizeof(blockScalar)==4 && sizeof(scalar)==8` and
also records the literal `sizeof(solveScalar)` with a reference to this entry.

**Impact.** Spec intent is met exactly; the native `simpleFoam` reference runs
with the same double-precision fields as `coupledFoam` (spec 3, "identical
precision"). No SPDP build is made.

## D-002 - Section 5.5 boundary-condition table is truncated in the spec

**Status:** decided.

**Finding.** The table was authored in Markdown with `|S_f|` inside cells; the
pipes split the cells. Lost content: the continuity-row entry of
fixedValue-U/zeroGradient-p, the wall diffusion coefficient, the p-p outlet
entry, the slip continuity entry and the totalPressure formula tail.

**Decision.** All BCs are made implicit in one uniform way, using the native
linearisation that every OpenFOAM patch field already provides:
`U_f = valueInternalCoeffs(w) . U_P + valueBoundaryCoeffs(w)` (componentwise)
and the same for `p`. Hence:

| Patch (U / p) | momentum rows | continuity row |
|---|---|---|
| fixedValue / zeroGradient | native `fvm::div` upwind + `fvm::laplacian` internal/boundary coeffs (wall: `a_P += nu_eff |S_f| deltaCoeff`, `b_P += nu_eff |S_f| deltaCoeff U_b`); pressure: `+S_f` on column 3 of the diag block (`p_f=p_P`) | `b^(3) -= U_b . S_f`; no p-p term (`D_f=0`, see D-013) |
| zeroGradient / fixedValue (outlet) | convection `a_P += phi_f`; pressure `b_P -= S_f p_b`; no diffusion | `+S_f^T` on U columns; p-p: `a^(3,3) += D_f |S_f| deltaCoeff`, `b^(3) += D_f |S_f| deltaCoeff p_b` (sign consistent with the `-laplacian` form, D-003) |
| slip / symmetry | native basicSymmetry coefficients (normal component removed) | `U_f . S_f = 0` -> no contribution |
| inletOutlet / totalPressure | native mixed coefficients (valueFraction from sign of phi_f); totalPressure `p_b = p0 - 1/2 |U_P|^2` lagged by one outer iteration | as above per face |
| cyclic / cyclicAMI | explicit-coupled: neighbour values lagged (source) | same |
| processor | implicit via block interface coefficients | same |

**Impact.** No BC is dropped. inletOutlet/totalPressure use native patch
evaluation, so any BC derived from fixedValue/zeroGradient/mixed works.

## D-003 - Section 5.3(d): dimension and sign of the p-p block

**Finding.** 5.3(d) writes the coefficient as `-D_f |S_f|^2 / |d_f . n_f|`;
5.3(e) uses `|S_f| / |d_f . n_f|`. Only the latter is dimensionally
consistent with the flux (m^3/s). The outlet cell of 5.5 reads
`a^(3,3) -= D_f ...`, which contradicts an M-matrix with negative
off-diagonals.

**Decision.** p-p off-diagonal `-D_f |S_f| deltaCoeff_f`
(`deltaCoeff_f = 1/|d_f . n_f|`, native `nonOrthDeltaCoeffs`), diagonal the
negative row sum plus the fixed-p boundary term `+D_f |S_f| deltaCoeff_f`.
Row 3 is the discrete `sum_f phi_f = 0`, so matrix and flux update (5.3e) are
built in the same face loop.

## D-004 - Test-doubleReduce reference value

**Finding.** 6.4 asks for "sum of 10^8 values of 1e-4 in float ... double result
must equal 1e4 to 1e-9 relative". `1e-4` is not representable in float:
`float(1e-4) = 9.99999974737875e-05`, so the exact sum of 10^8 such floats is
`9999.99974737875`, 2.5e-8 away from 1e4. The literal criterion fails for a
perfect accumulator.

**Decision.** The test checks the accumulator, not the representation: the
double-accumulated sum must equal `1e8 * double(float(1e-4))` to 1e-9 relative
(threshold unchanged). It also reports the difference to 1e4 and the float
result (expected wrong by >1 %, documenting the rule).

## D-005 - Section 12.3: "stop when both: (i) ... or (ii) ..."

**Decision.** Keyword `convergence.mode` (default `any`): the solver stops when
(i) `R_n < residualTol`, or, if a forceCoeffs function object is present,
when (ii) holds. `mode all` requires both. The benchmark harness evaluates (ii)
itself with the solvers' internal criteria disabled (spec 14), so the default
does not affect benchmark numbers.

## D-006 - MRF: `Omega x U`, not `2 Omega x U`

**Finding.** 5.6 says "Coriolis term 2 Omega x U ... relative frame formulation
as native MRF.DDt". Native MRF (`MRFZone::addCoriolis`) solves for the absolute
velocity and adds `Omega x U_abs`, which equals the relative-frame
`2 Omega x U_r + Omega x (Omega x r)`. The solver's `U` is absolute (native
fields and BCs, `MRF.correctBoundaryVelocity`), so `2 Omega x U` would be wrong
by a factor.

**Decision.** Implicit `V_P [Omega]_x` on the 3x3 momentum diagonal block,
i.e. exactly native `MRF.DDt` made implicit. Omega per cell is extracted from
the native `MRFZoneList::DDt` applied to the three unit vectors (public API
only; handles several zones and time-dependent Omega).

## D-007 - Hardware deviations

- WSL sees 94 GB RAM, spec assumes 128 GB. The 45 M-cell run (Phase 2) budget
  of 100-115 GB cannot be verified on this machine.
- SMT is on (32 logical / 16 physical cores). All MPI runs use
  `mpirun --bind-to core --map-by core`.
- The machine is shared with other long jobs; builds use `nice -n 19` and
  `WM_NCOMPPROCS=8`, which lengthens build time but not measured solver times.
  Benchmarks are only run when no other solver is running (checked and
  recorded in each benchmark JSON: `loadavg_before`).

## D-008 - Constants that the spec states but Section 11 does not list

Rule 0.5 forbids magic numbers; Section 11 is "complete" but several numbers
appear elsewhere. They become dictionary entries with the spec value as
default, under `coupled.guards` / `coupled.linear`:
`pivotGuard 1e-30` (6.3), `clampValue 1e30` (9.2), `maxRestarts 3` (6.2),
`maxCopAttempts 3` (6.3), `sentinel.UFactor 10`, `sentinel.pFactor 50`,
`sentinel.cflFactor 0.25` (9.3), `sentinel.maxReport 20`,
`remediation.warnFraction 0.01`, `remediation.warnInterval 100` (8.2),
`convergence.startupSwitchR 1e-2` (7.4), `nutMaxFactor` (9.2, listed).
The defaults are printed with the effective settings (11).

## D-009 - MRF has no gating test

**Finding.** Phase D implements implicit MRF, but no test case T0-T5 contains
an MRF zone, so the feature is not verified by any gate.

**Decision.** Add informational test `T-mrf` on the
`simpleFoam/mixerVessel2D` tutorial (MRF, SST): coupledFoam vs simpleFoam,
torque on the rotor within 2 %. Informational only (not a gate), because the
spec defines no threshold; the result is reported.

## D-010 - Build directory names

`~/OF/OpenFOAM-v2606-DP-Opt` and `~/OF/OpenFOAM-v2606-DP-Debug` instead of the
spec's `OpenFOAM-v2606-Opt/-Debug`: the precision is part of the name so a
DP and an SP/SPDP tree can never be confused.

## D-011 - Source file extensions

The spec's layout uses `.C/.H`; v2606 itself moved to `.cxx/.txx`. The project
keeps `.C/.H` as specified; wmake handles both.

## D-012 - Residual `b - A x` in double

The increment form (5.7) only helps if the right-hand side `b - A x` is exact.
It is evaluated in double from the double-precision assembly, before the
coefficients are narrowed to `blockScalar` (5.2 "conversion happens once per
outer iteration at store time"). The float matrix is only used for the
increment solve. This is what allows `R < 1e-8` (T0) with a float solver.

## D-013 - Rhie-Chow on boundaries

`D_f = 0` on boundary faces where p is not fixed (zeroGradient walls/inlets),
otherwise the "average gradient" term would create a spurious flux through
walls. On fixed-p faces `D_f = D_P` and the gradient term uses the cell
gradient. Coupled faces use the interpolated `D_f`.

## D-014 - Momentum diffusion assembled by the solver, not `divDevSigma`

`turbulence->divDevSigma(U)` would take the non-orthogonal correction limiter
from `fvSchemes`, which cannot vary per cell. To honour `nonOrthLimiter` (5.3a)
and `remediation.static.nonOrthLimiter` (8.1), the assembler builds
`-fvm::laplacian(nuEff, U)` with the uncorrected scheme plus an explicit
per-face-limited correction (native `limitedSnGrad` correction, the stronger
limit on faces touching a static cell), plus the explicit transpose term
`-fvc::div(nuEff dev2(T(grad U)))` exactly as native `linearViscousStress`.
The same limiter treatment is applied to the p-p Laplacian.

**Addendum (implementation):** the flux convention follows native
simpleFoam: `phi` is stored *relative* in MRF zones (`MRF.makeRelative`
after the flux update), because the turbulence model convects with the same
`phi` object. Spec 5.3(e) "makeAbsolute after the flux update" would make the
turbulence transport inconsistent in MRF zones.

## D-015 - blockILU0 is the block DILU form

"LDU-native incomplete block factorization with 4x4 block pivots" (6.3) is
implemented as the block analogue of native `DILUPreconditioner`: the
off-diagonal blocks of A are kept, only the diagonal blocks are modified,
`D*_u = D_u - L_f D*_l^-1 U_f`. For LDU (face-based) sparsity this is the
factorisation whose fill stays inside the LDU pattern; full ILU(0) would also
update off-diagonal blocks between two neighbours of a common cell, which
are not faces in LDU addressing and hence not stored.

## D-016 - Residual definitions used in the log

- `R` (7.1, CFL control) uses the L2 norm of 6.1: `||b - A x||_2 / normFactor`
  of the scaled system, divided by its value at iteration 1.
- `rU`, `rp` in the `CF|` line are L1-normalised like native OpenFOAM initial
  residuals: `sum|r_i| / sum(|Ax_i| + |b_i|)` over rows 0-2 and row 3. This
  keeps them directly comparable with the `Initial residual` of simpleFoam in
  the report plots.

## D-017 - Atomic write at file level

Spec 10 asks for fields to be written to `<time>.tmp` and renamed. A
directory rename conflicts with native Time bookkeeping: `uniform/time`
stores the directory name, and function objects keep state in the time
directory. Implemented instead: native field write, then `coupledState`
is written as `coupledState.tmp` and renamed (atomic `rename(2)`) per
processor directory. Because it is written last, its presence marks a
complete time directory, and restart reads only a directory that has it.
With the collated file handler the rename is skipped (logged); the case
templates use the uncollated handler.

## D-018 - "gradient limiter forced to cellLimited 1" in static cells

Implemented for the pressure gradient of the Rhie-Chow terms (native
`cellLimited Gauss linear 1` gradient selected per static cell). The
velocity gradient inside the high-order convection scheme comes from the
user's `div(phi,U)` scheme and cannot be overridden per cell without a
custom scheme. With the default `remediation.static.beta 0` static cells
use pure upwind convection, so that gradient does not enter there. If a user
sets `static.beta > 0`, the unlimited gradient is used in those cells. This
is logged as a limitation.

## D-019 - Outflow through a fixedValue velocity face

The 5.5 table (as far as legible) gives `a_P += max(phi_f, 0)` for fixedValue
faces. The implementation uses the native `fvm::div` boundary treatment,
which takes the boundary value for both flow directions (`b -= phi_f U_b`).
The two differ only if the flux through a fixedValue face points against the
prescribed velocity, which does not occur in the test cases.

## D-020 - A diverged linear solve is a failed step

Observed on T0 at CFL 500 (2026-09-21): near the float floor block-BiCGStab
can diverge (final linear residual 9e4 after 200 iterations). The spec's
rule "after maxCflCuts accept the step with omegaMin" then applies 10 % of a
garbage increment and destroys a converged state. Two changes:
1. blockBiCGStab returns the iterate with the smallest residual if it does
   not converge (Krylov output never worse than the initial guess).
2. A solve that does not reduce the residual, or yields non-finite values,
   is treated like `omega < omegaMin`: CFL cut by kappa and the step is
   repeated. If the cuts are exhausted, the step is skipped (no field
   update) and CFL is cut by the sentinel factor. The line-search rule of 7.2
   is unchanged for successful solves.

## D-021 - Strong pressure reference in closed domains

For a closed domain (`p.needReference()`), the native weak reference
(doubling the diagonal of the reference cell) leaves a near-null space. After
the right-hand side is rounded to float, its compatibility condition is
violated slightly and the float Krylov solve drifts along the
constant-pressure mode. Observed on T0: dp-limited line search, then
divergence. Instead, the continuity row of the reference cell is replaced by
`d (p_ref - pRefValue) = 0`. This is exact, not an approximation: in a
closed incompressible domain the continuity equations sum to zero
identically, so one of them is redundant. Result on T0: stable convergence to
R = 7e-8.

## D-022 - Linear tolerance vs. the outer residual target

Spec 6.1 sets the linear `tolerance 1e-8` (absolute, on ||r||/normFactor);
T0 requires the outer `R < 1e-8`, where R is normalised by the iteration-1
value (||r_1||/normFactor = 0.088 on T0). Once ||r||/normFactor < 1e-8 the
linear solver does zero iterations, the increment is zero and the outer
iteration stalls (observed: R frozen at 7.0e-8). The T0 case therefore uses
`tolerance 1e-12` (a tightening, allowed by rule 0.2). The default stays at
1e-8; cases with a tighter outer target set the linear tolerance
accordingly.

## D-023 - Test-blockGAMG solve tolerance 1e-9

Spec 6.4: "solve to 1e-8; iteration count <= 20; identical solution 1 vs 4
ranks to 1e-5". With the 6.1 norm (||r||_2 / normFactor, normFactor an L1
sum over all rows) 1e-8 is only a ~5e-6 relative L2 reduction, and two
solutions stopped there differ by 1.3e-5 (measured, cells matched by
centre), so the cross-rank criterion fails although both solves are correct.
Measured on the T0 mesh:

| tolerance | iterations np1/np4 | ||x1-x4||/||x1|| |
|---|---|---|
| 1e-8  | 9 / 10  | 1.32e-5 |
| 1e-9  | 11 / 13 | 9.3e-7  |
| 1e-10 | 13 / 14 | 6.0e-7  |
| 1e-11 | 15 / 16 | 6.0e-7 (float floor) |

The test uses 1e-9: a tightening (rule 0.2) that meets both remaining
criteria. The cross-rank metric is the relative L2 difference of the whole
solution vector; per-component integrals are recorded too but not used,
because the pressure-component integral nearly cancels (~2e-6 against ~8e-5)
and its relative difference is not a measure of solution identity.

## D-024 - "Outer iterations to R < 1e-5" for simpleFoam (T2)

simpleFoam has no combined residual R. For the T2 iteration-count criterion
the native iteration count is the first iteration at which all initial
residuals (p, Ux, Uy, k, omega) are below 1e-5; for coupledFoam it is the
first iteration with R < 1e-5. The two normalisations differ; the metric is
reported as defined here, and the report additionally compares both solvers
with the identical force/pressure-drop window criterion of 12.3(ii) and in
wall-clock and CPU time, which are the decisive quantities for the user.

## D-025 - Benchmark configurations A and B, time metrics

- A "simpleFoam with tutorial settings": the case's own native solver
  settings and relaxation factors as shipped with the tutorial (for
  pitzDaily and backwardFacingStep2D the tutorial itself uses SIMPLEC,
  consistent yes).
- B "SIMPLEC, p 1.0 / U 0.9 / k,omega 0.9" as specified.
- Stop: the solvers' own criteria are disabled (fixed iteration budget);
  the harness finds the convergence iteration with the window criterion.
- Time to convergence: every rank runs under `/usr/bin/time -v`; wall time
  and CPU time (user + sys summed over all ranks, reported in CPU-hours) of
  the whole run are scaled by the solver's own time progression at the
  convergence iteration (coupledFoam `tWall`, simpleFoam `ClockTime`).
  Speed-up is reported both as wall-time and as CPU-hour ratio (user
  request 2026-09-21: the fastest solver, not the one with the fewest
  iterations, is the goal).

## D-026 - Anderson acceleration: Walker & Ni Type II update (B5)

The short formula in B5, x_(k+1) = x_k + beta (w_k - sum alpha_j dw_j) with
x_k taken after the update, omits the state-difference term and, read
literally, adds w_k twice. Implemented is the Type II form of Walker & Ni
(2011): with x_k the state the step started from and f_k = omega dx the
accepted increment,
    gamma = argmin || S (f_k - sum_j gamma_j dF_j) ||_2
    x_(k+1) = x_k + beta f_k - sum_j gamma_j (dX_j + beta dF_j),
S scaling velocity by Uref and pressure by pref (dimensionless least
squares). The least-squares problem is solved with an updated QR
factorisation (two Gram-Schmidt passes, Givens down-dating), all dot
products in double and global. Memory is 2m+2 double 4-vectors per cell
(m = 4: 320 B/cell, 14.4 GB at 45 M cells, slightly above the 11.5 GB of the
B7 table). The Anderson history is not part of the restart state (a
restart starts with an empty history). Flushes: CFL cut, skipped step,
omegaMin marking, dynamic-set size change, rollback, sentinel rejection.

## D-027 - Zonal factors act only through dt and beta (B6)

Zonal factors act only through the local pseudo-time step (cflFactor) and
the convection blending factor (beta). A smaller local dt adds a larger
diagonal term V/dt, i.e. implicit local under-relaxation; the Newton-like
update stays consistent. Explicit per-patch or per-zone relaxation factors
on the increment are deliberately not implemented: they would change the
direction of the increment that the line search judges. Processor
boundaries need no relaxation (block interfaces are fully implicit).
Patch-distance layers are computed once by breadth-first growth from the
patch faces' cells, across processor faces; processor patches are never
start patches. Several matching entries multiply; the zonal factors apply
after the static/dynamic rules (dynamic cells keep beta = 0). Input checks:
cflFactor > 0, 0 <= beta <= 1, nLayers >= 1.

## D-028 - The block Gauss-Seidel smoother is not adequate (measured)

T0 Re 1000 (128x128, linearUpwind deferred correction, 300 outer
iterations, 1 rank, before amendment B):

| linear solver | iterations to R < 1e-8 | CFL cuts | sum linear its | t_solve |
|---|---|---|---|---|
| BiCGStab + blockGaussSeidel | not reached (R 3.4e-4) | 29 | 13757 | 950 s |
| GMRES + blockGaussSeidel | not reached (R 6.9e-2) | 0 | 58716 | 811 s |
| BiCGStab + blockILU0 | 59 | 0 | 406 | 26 s |
| GMRES + blockILU0 | 59 | 0 | 374 | 7 s |
| BiCGStab + GS, pure upwind (beta 0) | 186 | 14 | 9106 | 509 s |

The outer-iteration oscillation at Re 1000 was first attributed to the
explicit deferred correction; the ILU0 runs disprove this: with the block
DILU smoother the same discretisation converges in 59 iterations at CFL
500 without a single cut. Poor increments from the block Gauss-Seidel
smoothed solve (the 4x4 point smoother is weak on the saddle-point
structure of the coupled system) caused the oscillation. The spec's formal
ILU0 decision (Gate B2, >= 15 % wall-clock gain on T2) is taken on T2; the
T0 data above are the supporting evidence. A residual-growth rollback that
was considered for the supposed deferred-correction instability is not
implemented (not needed).

## D-029 - Linear-solve failure: B4 replaces D-020 item 2

With amendment B4, a linear solve that does not reach eta*||r0|| (or the
absolute floor) within maxIter is a failure: CFL is cut by kappa and the
step repeated (counting toward maxCflCuts); maxLinFails consecutive
failures abort through the 9.3 diagnostic path (last-valid fields written,
FatalError). D-020 item 1 (best-iterate return) stays. The "skip the step"
fallback of D-020 remains only for the case maxCflCuts < maxLinFails.

## D-030 - Processor agglomeration (6.3.4): keywords, mechanism, level 0

`blockGAMG.processorAgglomerator`: `masterCoarsest` (default) applies rule
6.3.4 - T(N) = nRanks, max(1, nRanks/procAgglomDivisor) for
N < procAgglomCellsPerRank*nRanks, 1 for N < procAgglomCellsPerRank -
through a native GAMGProcAgglomeration subclass (runtime name
`blockGAMGRule`, header blockGAMGProcAgglomeration.H); `none` disables
processor agglomeration; `nativeMasterCoarsest` selects the native class;
any other word is passed verbatim to the native selector (`manual`,
`procFaces`, `eager`). Parameters `procAgglomCellsPerRank` (5000) and
`procAgglomDivisor` (4) from coupledDefaults. Groups are contiguous,
balanced rank blocks (master = lowest rank), like the native
masterCoarsest. The coarsest level is gathered by blockGAMG itself
(lduPrimitiveMesh::gather; the native mechanism cannot gather the coarsest
level alone). Level 0 is never agglomerated - a native restriction - so
ranks-per-level reports nRanks at level 0 even where the formula asks for
less; the unit-test rule check excepts level 0.

Verification (run/unit_cavity, 16 k cells): serial bit-identical to the
pre-change build; 4 ranks on/off relDiff <= 2.0e-6 across V/F/W cycles,
two-stage (4->2->1) and multi-stage (4 4 2 1 1) gathers and the dense-LU
coarsest solve; `nativeMasterCoarsest` pass-through works. Known issue
found during this verification: the K cycle stalls just above tolerance
on unit_cavity (1.8e-8 vs 1e-9) in the pre-change baseline too -
pre-existing, tracked in HANDOFF 4.1, unrelated to 6.3.4.

## D-031 - User directives of 2026-09-21 (resources, meshes, deliverables)

- Heavy tests use at most 10 cores: CF_HEAVY_NP default is now 10 (was
  16), with MPI core binding on a free machine. The earlier arrangement
  next to the F1 job (CF_MPI_CPUSET=10-15, CF_HEAVY_NP=6) is obsolete
  once that job is stopped; CF_MPI_CPUSET remains available.
- T-scaling therefore runs ranks 1, 2, 4, 8, 10 by default instead of the
  spec's 1, 2, 4, 8, 16 (recorded as specRanks=false in the result; the
  0.8 relative-efficiency pass threshold is unchanged).
- T4 mesh variant b is retargeted from 3-5 M to 1-2 M cells (surface
  level (6 6), features 7, refinementBox 5): the user wants a presentable
  motorBike case, reusable as demonstration material for their
  post-processing tool. Not a threshold change; all pass criteria stay.
- Every snappyHexMesh mesh is built once (run/T4a_mesh, run/T4b_mesh,
  run/T5_mesh) and reused by all tests and benchmarks - this was already
  the design (test_T4_motorBike.mesh_dir) and is now an explicit user
  requirement.
- The paper is delivered as a compiled PDF (LaTeX installed 2026-09-21:
  pdflatex, latexmk, bison present), with intermediate PDF builds during
  the work so the user can inspect drafts.

## D-032 - Two paper versions: professional and tutorial (user, 2026-09-21)

The user wants TWO final PDFs built from the same results/ data and the
same generated figures/tables (bench/make_report.py):

1. `report/paper/paper.pdf` - the professional paper as before.
2. `report/paper/paper_tutorial.pdf` - a tutorial version for readers
   without a CFD background ("for dummies", the user\x27s words): same
   structure and numbers, but every concept explained in plain language
   (what a coupled vs segregated solver is, why pressure-velocity
   coupling is hard, what multigrid/preconditioning/CFL/pseudo-time
   stepping do, what the test cases show and why the reader should
   care), with intuition boxes around the equations rather than more
   equations. The user wants to learn from it and understand what was
   done and why.

Both share numbers.tex, figures/ and tables/ so neither can drift from
the measured data; the Makefile builds both. Intermediate PDF builds of
both are delivered to the user during the work (D-031).

## D-033 - Benchmark accounting and initialisation (review findings, 2026-09-21)

- D-025 addendum. The cost of coupledFoam includes its potentialFoam
  initialisation (Allrun, `coupled.potentialInit yes`). potentialFoam
  runs under the rank wrapper, and its wall and CPU time are added to the
  totals; both parts are recorded (potentialFoamWallSeconds /
  potentialFoamCpuHours, solver* / pre*). Time to convergence is
  t_pre + f*t_solver, and the same for CPU, where f is the solvers

## D-034 - Restart state per processor, ASCII (review fix)

`coupledState` is a per-processor `localIOdictionary`, always written as
ASCII with 17 significant digits. It holds:
- CFL/PTC state
- Eisenstat-Walker state
- the dynamic set
- the gamgAutoTune state (sweeps, cycle type, rho window, hysteresis,
  failure counters, latches), re-applied to the GAMG on restart
- the last Cd/Cl window of the convergence monitor

`coupledD` and `coupledQ` are written next to it. The Anderson history
is not saved (D-026).

Before this fix:
- The state was a global IOdictionary, so in parallel every rank read
  rank 0's state. On T0 np4 the first residual after a restart jumped
  from 5.6e-5 to 1.3e-2.
- A binary dictionary cannot read back an empty list `0 ( )` (restart
  crash: "dynamicSet has 2 excess tokens").

Verified: the restart round trip is bit-exact on T0 np1 (30 + 26 = 56)
and np4 (30 + 46 = 76, same final R as the uninterrupted run on the same
decomposition). Old binary coupledState files cannot be read.

## D-035 - Conservative flux; Rhie-Chow boundary weighting (review fixes)

- After the solve, the flux update reuses the explicit Rhie-Chow q_f of
  the solved assembly (A7). Sum phi per cell then equals that continuity
  row's linear residual, as in native `phiHbyA - pEqn.flux()`. q_f is
  stored as `coupledQ` for the spec-10 restart check.
- On physical patches the explicit q_b is weighted by
  (1 - valueInternalCoeffs(p)), the same weight as the implicit p-p term
  (A3). This matters for mixed-type BCs (freestreamPressure,
  inletOutlet): before, their zero-gradient faces carried a spurious flux
  that D-013 was meant to prevent.
- The non-orthogonal part of q uses the cellLimited gradient in static
  cells (A8, extends D-018).
- Processor faces decide the non-orthogonal limiter by the same rule as
  internal faces: static if either cell is static. The static flag is
  exchanged across processor patches (A2).
- Wedge patches go through the general BC path instead of being treated
  as zero-flux (A4).
- Row 3 is scaled by s_p in double before narrowing to float (A5).
- normFactor excludes the pressure-reference row and reduces on
  mesh.comm() (A6).
- processorCyclic interfaces are not block-coupled (A9).
- A closed domain whose boundary fluxes do not balance gets a one-time
  Warning; tolerance refFluxBalanceTol = 1e-8, the native adjustPhi value
  (A10).
- The nonOrthCorrection constructor does one global reduction instead of
  one per patch; per-patch reductions caused the MPI_ERR_TRUNCATE (f025546).

## D-036 - Line-search beta (spec 7.2)

beta >= 1 multiplies CFL after a full step with no cuts, capped at CFLmax
and not applied during a hold. The default of 1 is spec 7.2's "no
additional boost". Before, beta was read but unused.

## D-037 - autoTune demotion reachable (6.3.5); Anderson flushes; aborts

- A low window at nPostSweeps == 1 with rho < 0.2 after a controller
  promotion now demotes one step. Hysteresis is a counter
  (tuneConsecutiveWindows).
- The Anderson history is flushed on:
  - a CFL cut or skipped step
  - the line-search minimum
  - a rollback
  - a change in dynamic-set membership (version counter, not set size)
  - the end of start-up (beta 0 -> 1)
  - any autoTune change
- Skipped steps no longer update Eisenstat-Walker, PTC or the sentinel.
  The rho window takes only accepted, converged solves. The convergence
  monitor does not record a repeated (Cd, Cl) pair.
- Both abort paths write fields and remediation sets only into
  `<iter-1>_lastValid`, named after the iteration whose fields it holds.
  No more incomplete time directories.
- Anderson is disabled with a Warning above andersonMaxCells = 35 M
  cells (B7). maxLinFails < 1 is a FatalIOError.
- blockGMRES gets the same mixed-precision refinement as blockFGMRES
  (double iterate, true double residual).
- potentialFoam `-writep` needs `div(div(phi,U))`; it was added to the
  fvSchemes of T1-T5, where every coupledFoam run failed at potentialFoam
  before.

## D-038 - T1 study: linear-solver defaults and T1 k/omega scheme (2026-09-22)

Study: 17 serial runs on T1 (pitzDaily, kOmegaSST, 12225 cells). The
logs are in run/exp_T1_*; the summary is run/exp_T1_results.txt.

Findings:
1. With the tutorial's second-order k/omega convection (limitedLinear 1)
   T1 has no convergent steady state for EITHER solver:
   - simpleFoam's p residual stalls at 6e-4 and dp wanders by +-0.07 %.
   - linearUpwind k/omega stalls at 1e-3.
   - With upwind k/omega simpleFoam converges to 1e-8 in 769 iterations
     (72 s).
2. FGMRES restart 10 caps CFL: GMRES(10) stagnates, the linear solve
   fails at CFL ~50-70, and each failure cuts CFL and holds it. With
   restart 30, CFL reaches 500 without cuts. Block Gauss-Seidel keeps CFL
   at 1-2 (confirms D-028).
3. An absolute linear tolerance of 1e-8 cannot reach R < 1e-6 when the
   normaliser R1 is 6e-3. The linear tolerance must stay below
   R1*residualTol, as for T0 (D-022).
4. Eisenstat-Walker is not the cause. Fixed relTol 0.1 or 0.01 was no
   better, and with frozen turbulence tighter solves made R grow.
5. The block-GAMG preconditioner is weak on graded, high-aspect-ratio
   meshes. rho (the first application) is often > 1, and a 0.5 reduction
   needs 20-60 FGMRES iterations. This leaves coupledFoam about 5x SLOWER
   than simpleFoam in wall time on T1 (converging configuration: 483 its,
   364 s; simpleFoam 72 s). Open: FABLE_REVIEW.md item 2.

Decisions:
- T1-T5 defaults: smoother blockILU0, FGMRES restart 30 (B7: restart 6
  above 40 M cells stays) and linear tolerance 1e-10. T0 is unchanged
  (validated).
- T1 uses `bounded Gauss upwind` for div(phi,k) and div(phi,omega) for
  BOTH solvers, and relaxes k/omega by 0.95 in coupledFoam.
  - Spec 14 still holds: both solvers use identical schemes.
  - The pass criteria are unchanged; they compare against simpleFoam on
    the same case files.
  - The upwind k/omega answer differs from the limitedLinear one by
    0.95 % in dp (stalled limitedLinear simpleFoam as the reference).
  - The paper reports this deviation from the tutorial.

## D-039 - Preconditioner study: matrix-weighted agglomeration and one finest ILU0 sweep (2026-09-22)

Study on branch `precond-research` (worktree /home/jonas/cf_precond, own
install /home/jonas/cf_precond_platform), answering FABLE_REVIEW item 2.
Method: `coupled.dumpLinearSystem (iters)` writes the linear system of an
outer iteration; `Test-blockSystem` replays it with any number of solver
variants (iterations, rho, the scale-free rhoOpt = min_a ||r - a A M^-1 r||
/||r||, set-up and apply time). Systems from T1 iterations 30, 150, 300 and
450 (CFL 98, 58, 500, 500).

Findings (T1, np 1):
1. The diagnosis "rho > 1" was partly a measurement artefact: rho as
   defined in 6.3.5 is unscaled, so a preconditioner that is right in
   direction but wrong in magnitude scores > 1. The scale-free rhoOpt of
   the committed configuration was nevertheless 0.99-0.9999, i.e. one
   cycle removed almost nothing of the residual.
2. The finest-level smoothing was the main defect. blockILU0 as a
   Richardson smoother (x += M^-1 (b - A x)) DIVERGES on this saddle-point
   system: the rows are far from diagonally dominant (max |A_kc| of the
   neighbours over |D_kk|: 9.5 in the momentum rows through the pressure
   column, 8.7 in the continuity row through the velocity columns). With
   nFinestSweeps 3 instead of 2 the solve needs 21 instead of 16
   iterations; with 1 sweep it needs 1. Block Gauss-Seidel diverges
   outright here (rho 1e7).
3. Geometric agglomeration (faceAreaPair) on the graded mesh is the second
   defect. Pair agglomeration on block-matrix weights
   w_f = |a_uu,f|/sqrt(a_uu,P a_uu,N) + |a_pp,f|/sqrt(a_pp,P a_pp,N)
   ("combined") gives rhoOpt 0.33-0.45 instead of 0.99 at CFL 500.
   Momentum-only or pressure-only weights are worse than the sum.
4. autoTune makes it worse: it reads the unscaled rho (~1.1), calls it bad
   and raises nPostSweeps to 4, which costs 40 % wall time for nothing.
5. A SIMPLE-type block preconditioner (blockSimple, S = schurScale*C, one
   blockGAMG hierarchy on the decoupled [A 0; 0 S]) is the best
   preconditioner per application (rhoOpt 0.59-0.68, 3 FGMRES iterations
   in sequential mode) but costs two cycles, and the outer iteration
   DIVERGED on T1 (R 1.8 after 600 iterations): its velocity increment is
   not accurate enough for the line search. Kept as a selectable option,
   not a default.
6. Equilibration (symmetric diagonal scaling of the block system) does not
   help; CFLmax 2000/10000 does not reduce the outer iteration count
   (529/540 vs 527); etaMax 0.1 reduces outer iterations 527 -> 432 but
   costs more wall time.

T1 wall clock (np 1, same machine load, R < 1e-6):

| configuration | outer its | wall s | mean lin its | mean rho |
|---|---|---|---|---|
| committed (D-038) | 518 | 175 | 25.4 | 12.7 |
| + combined weights | 521 | 77 | 7.8 | 4.3 |
| + nFinestSweeps 1 | 491 | 44 | 4.0 | 1.1 |
| + nCellsInCoarsestLevel 20 | 527 | 41 | 4.0 | 1.4 |
| + autoTune no + reagglomerateInterval 50 | 495 | 33 | 3.3 | 1.0 |
| simpleFoam, same schemes, to 1e-8 | 769 | 37 | - | - |

Decisions:
- New selectable keywords in the blockGAMG dictionary, all defaulting to
  the previous behaviour: `agglomerationWeights`
  (geometric | momentum | pressure | combined), `reagglomerateInterval`,
  `scaleCorrection` (none | finest | all). New preconditioners
  `blockSmoother` (single level, for diagnosis) and `blockSimple`.
- T1 (and the T2-T5 templates that follow it) use agglomerationWeights
  combined, reagglomerateInterval 50, nFinestSweeps 1,
  nCellsInCoarsestLevel 20 and autoTune no. T0 is unchanged: it uses the
  blockGaussSeidel smoother, which needs its two finest sweeps
  (nFinestSweeps 1 there: 186 instead of 56 outer iterations).
- The library defaults stay as they are, because the evidence is from the
  ILU0 smoother on graded meshes only.
- Open: rho of 6.3.5 should be measured scale-free (rhoOpt) before
  autoTune can be trusted; the startup phase (iterations 1-60, upwind,
  CFL below 100) still needs 7-10 linear iterations per solve.
## D-040 - nut safety cap raised to 1e8*nu (T3 root cause, 2026-09-22)

The applyBounds nut cap (nutMaxFactor, spec 9.2) was 1e5*nu. On T3
airFoil2D (kOmegaSST, Re_c ~ 9e7) the correct solution has nut/nu up to
3.6e5 in the wake (measured on the simpleFoam reference), so the cap
clipped the wake eddy viscosity to 28 % of its physical value on every
outer iteration; the under-diffused shear layer went numerically
unsteady, R limit-cycled at ~1e-2 and the forces oscillated with sign
flips (the reported Cd -0.015 / Cl -0.46 were one sample of that
oscillation, not a sign error; the GEKO run on the identical case capped
almost no cells and was fine). nutMaxFactor is now 1e8: a divergence
guard that never binds on a physically correct RANS solution. The
per-case override coupled/bounds/nutMaxFactor is unchanged.

Verified: with the cap at 1e8, nNutCapped stays 0 on T3 (run
exp_T3_nutcap). NOT resolved by this alone: T3-SST still limit-cycles
(R 1e-3..1e-2) - tracked in FABLE_REVIEW (suspects: PTC aggressiveness
vs a case simpleFoam needs 20000 relaxed iterations for; k bounding
firing every iteration; freestream mixed-BC switching).

## D-041 - Hot loops are hand-written; expression templates out of scope (B11 6.5.7)

The block kernels (Amul, residual, smoothers, restriction/prolongation,
all Krylov vector operations) are hand-written loops over flat
blockScalar* arrays with __restrict__, `#pragma omp simd` (-fopenmp-simd,
no OpenMP runtime) and 64-byte aligned storage (alignedList), plus fused
kernels (axpy_dot, update_residual_norm, fused modified Gram-Schmidt)
that read each operand once (blockKernels.H). An expression-template
layer over block fields is explicitly NOT introduced: the dominant
kernels are a handful of already single-pass, bandwidth-bound loops
(gate phase_A: Amul 85.6 %, axpy_dot 141.5 % of the STREAM triad on a
5 M-cell system, results/gates/phase_A.json), OpenFOAM field algebra is
under 30 % of iteration time and lives in core, and a template layer
would hide the precision crossings D-001/D-012 require to be explicit.
The 4x4 matvec uses an explicit 4-wide vector product with an
in-register transpose that preserves the scalar summation order
(bit-identical); SIMD reductions reorder double partial sums, which is
accepted (D-012 keeps accumulation in double; unit tests hold at their
tolerances) - measured effect: T0 Re100 56 -> 57 outer iterations, i.e.
inside the trajectory scatter of FABLE_REVIEW item 1, for ~37 % less
solve time per linear iteration.

## D-042 - Averaged force comparison for wake cases (user decision, 2026-09-22)

T4a, T4b and T5 have physically oscillating wakes: the T4b
simpleFoam reference oscillates Cd by +-1 % and Cl by +-8 %
indefinitely (run/ref_T4b_np10), so criterion 12.3(ii) (0.2 %
min/max window) is unsatisfiable for any steady solver on these
cases, coupledFoam and simpleFoam alike. The user approved
averaging and coarser tolerances for exactly these cases ("bei den
motorbike und ahmed body cases muessen wir wahrscheinlich averagen
und das vergleichen und vllt auch groebere toleranzen nehmen wegen
dem nachlauf problem. das ist okay", 2026-09-22). New evaluation,
oscillatory cases only: window W = max(500, n/4) (<= n/2);
stationarity = the two half-window means differ by <=
max(0.5 % |mean|, 0.002) for Cd and Cl; comparison = window means,
Cd within max(2 % relative, 0.002 absolute), Cl within
max(2 % relative, 0.01 absolute) of the reference. Time to
convergence uses the first stationary window. T0-T3 criteria are
unchanged. This is a user-approved criterion change, not an agent
loosening.

## D-043 - Robust linear-solver defaults: damped ILU0 V cycle (TASK 1, 2026-09-22)

Problem: the D-039 settings (K cycle, combined weights, nFinestSweeps 1,
nCellsInCoarsestLevel 20) were tuned on T1 only. They caused B4 aborts on
T2, T3 and T4a. On the serial T4a mesh they are not even admissible: the
K-cycle coarsening-ratio rule r_l >= 3 (6.3.1) rejects every mergeLevels
2..4 with coarsest 20 (last-level ratio 2.09 -> FatalError). The
conservative set (V / geometric / 2 sweeps / coarsest 200) also degraded
on T4b at CFL ~30 (FABLE_REVIEW 4a).

Method (OPUS_TASKS TASK 1a-1c, applied mechanically):
- Linear systems dumped from:
  - T1: iterations 300 and 450 (CFL 500).
  - T2: iterations 400 and 800.
  - T3-SST: iterations 30 and 100. The run aborted via B4 at 125, so
    the last two existing dumps were used.
  - T4a serial: iterations 40 and 50. This dump run used the V cycle
    because the D-039 K/c20 set is inadmissible there, so these dumps
    sit at CFL 17/21, not 500.
- 12 variants solved by Test-blockSystem with eta = 0.1.
- PASS = converged to a 0.1 reduction within 30 FGMRES iterations AND
  rhoOpt <= 0.9.
- Table cells: its/rhoOpt; x = fail. Logs: run/robust_*.log. JSON:
  results/exploratory/task1_grid.json.

| variant | T1-300 | T1-450 | T2-400 | T2-800 | T3-30 | T3-100 | T4a-40 | T4a-50 | sum t1 [s] | all pass |
|---|---|---|---|---|---|---|---|---|---|---|
| V_geo_f2_c200 | 200/1.00 x | 200/1.00 x | 3/0.14 | 400/1.00 x | 55/0.93 x | 23/0.97 x | 12/0.46 | 5/0.46 | 20.92 | no |
| V_comb_f2_c200 | 11/1.00 x | 10/0.96 x | 2/0.17 | 400/1.00 x | 20/0.96 x | 30/0.87 | 8/0.41 | 5/0.47 | 15.31 | no |
| V_comb_f1_c200 | 7/0.85 | 6/0.58 | 3/0.34 | 15/0.95 x | 58/0.80 x | 23/0.67 | 15/0.70 | 54/0.75 x | 18.00 | no |
| V_comb_f1_c20 | 7/0.83 | 6/0.56 | 7/0.35 | 11/0.96 x | 58/0.80 x | 23/0.67 | 46/0.70 x | 400/0.77 x | 92.13 | no |
| V_geo_f1_c200 | 200/0.99 x | 200/0.96 x | 3/0.28 | 17/0.93 x | 400/0.70 x | 400/0.79 x | 22/0.71 | 27/0.71 | 24.22 | no |
| K_comb_f1_c20 | 8/0.42 | 6/0.33 | 18/0.43 | 12/0.28 | 62/0.80 x | 26/0.67 | FATAL | FATAL | 1.14 | no |
| K_comb_f1_c200 | 12/0.46 | 9/0.46 | 3/0.33 | 25/0.29 | 400/0.92 x | 400/0.94 x | 40/0.73 x | 74/0.77 x | 33.47 | no |
| K_geo_f2_c200 | 144/1.00 x | 119/0.99 x | 2/0.13 | 400/0.95 x | 400/0.98 x | 400/0.99 x | 12/0.43 | 7/0.47 | 28.82 | no |
| K_comb_f2_c200 | 26/0.98 x | 20/0.93 x | 2/0.14 | 400/1.00 x | 400/0.99 x | 27/0.99 x | 11/0.44 | 7/0.49 | 26.52 | no |
| GS_control | 200/1.00 x | 200/1.00 x | 400/1.00 x | 400/1.00 x | 400/1.00 x | 400/1.00 x | 400/0.84 x | 400/0.81 x | 180.52 | no |
| V_comb_f2_c200_r07 | 4/0.59 | 2/0.26 | 2/0.19 | 7/0.39 | 3/0.60 | 3/0.65 | 7/0.43 | 4/0.45 | 6.60 | YES |
| V_comb_f2_c200_r05 | 3/0.32 | 3/0.26 | 2/0.29 | 2/0.28 | 3/0.35 | 3/0.34 | 10/0.59 | 6/0.62 | 7.09 | YES |

No undamped variant passes all dumps: the T3 start-up dump at CFL 40
defeats all of them. TASK 1d was therefore executed:
- New keyword blockGAMG.smootherRelaxation for blockILU0 and
  blockGaussSeidel: x += relax*(smoothed - x) per sweep.
- Default coupledDefaults::smootherRelaxation = 1, the bit-identical
  undamped path; values are validated to lie in (0, 1].

Both damped variants pass every dump. By the selection rule (smallest
summed wall time) the winner is V_comb_f2_c200_r07.

Defaults applied to T1-T5 (T0 untouched):
- Changed: cycleType V, agglomerationWeights combined, nFinestSweeps 2,
  nCellsInCoarsestLevel 200, smoother blockILU0, smootherRelaxation 0.7.
- Unchanged: autoTune no, reagglomerateInterval 50, restart 30, linear
  tolerance, maxIter 400.

Interpretation:
- The undamped point-block ILU0 over-corrects on the non-dominant
  saddle-point rows (run/robust_*_dominance.log):
  - continuity row: offU/d up to 30-34 on T2/T3
  - momentum row: offP/d up to 15 on T1
- Damping restores smoothing without the second-sweep divergence that
  D-039 observed.
- The V cycle avoids the K-cycle ratio-rule failure on snappy meshes.

## D-044 - Dynamic-set thresholds cU 2 -> 4, cp 5 -> 15 (user, 2026-09-22)

The user (F1 aerodynamics) pointed out that the spec 8.2 defaults mark
physically correct cells.

Definitions (lineSearch::setReference):
- U_ref = max |U| over the field and the non-coupled boundaries, which is
  essentially the free stream.
- p_ref = 0.5 U_ref^2.

What the old defaults mean:
- cU = 2: |U| > 2 U_inf. Under wings with a strong suction peak the
  local speed reaches 2-3 U_inf.
- cp = 5: Cp < -5. The user has seen suction peaks of Cp -7 to -8.

Consequence: exactly the cells that produce the downforce were put into
the dynamic set. There they get:
- beta = 0 (first-order upwind)
- a local CFL factor of 0.1
- increment clipping
This is a silent loss of accuracy where it matters most.

New defaults, consistent through Bernoulli (|U| = cU U_inf gives
Cp = 1 - cU^2):
- cU = 4, which covers 2-3 U_inf with margin.
- cp = cU^2 - 1 = 15, which covers Cp down to -8 with margin.

These stay well below the sentinel (divergence) limits of 10 U_ref and
50 p_ref (sentinelUFactor, sentinelPFactor), so real blow-ups are still
caught first by the dynamic set and then by the rollback. cSpike (0.5,
the jump to the neighbour mean) is unchanged.

Applied to coupledDefaults.H and to the case templates T1-T5. T0 keeps
its validated explicit values 2 / 5 (lid-driven cavity, |U| <= 1, where
the thresholds never bind).

### D-043 addendum - smootherRelaxation 0.5 instead of 0.7 (same day)

Acceptance check 1f on T2 (pytest, D-043 settings): B4 abort at outer
iteration 296 (CFL 500, R 3.7e-4). The T2 grid dumps above came from the
D-039 trajectory and reached only CFL 9 / 47, so the CFL-500 regime was
not covered. New dumps run/dump_T2b iter 285 / 293 (CFL 500): no
variant passes the 1c rule there (preconditioner amplifies the residual,
||Az-r||^2/||r||^2 ~ 408 for r05; worst cells distributed over the wake,
top-10 share 12 %). Note: dump_T2b was created from the D-043 template,
so every variant without an explicit smootherRelaxation inherited 0.7 on
those two dumps (V_comb_f2_c200 == r07 there).
Decisive full-run probe (T2, 1000 outer its, serial, run/exp_T2_*):

| variant | wall | end R | B4 |
|---|---|---|---|
| D-043 winner (relax 0.7) | - | - | abort at 296 |
| relax 0.7 + CFLmax 100 | 418 s | 2.5e-5 | none |
| relax 0.5 | 259 s | 8.9e-6 | none |
| damped GS 0.7, V/geometric/f2/c200 | 157 s | 2.9e-5 | none |

relax 0.5 passed all eight original dumps as well (sum 7.09 s vs 6.60 s
for 0.7). Chosen default for T1-T5: smootherRelaxation 0.5 (lowest R,
no abort). The CFL-500 preconditioner weakness on T2 (and the T4b
degradation at CFL ~30 with the undamped set) stays open in
FABLE_REVIEW item 4.

### D-042 addendum - larger window, coarser stationarity, mean-field delta comparison (user, 2026-09-22)

User decision (2026-09-22, translated from German): "T4b needs a coarser
criterion. It will probably not converge even with more iterations
because the flow is simply that unsteady. So maybe we just need a larger
averaging window and different tolerances. And for the comparison maybe
also look at the flow field as a delta and compare it, not just numbers."

Evidence (cached references, read-only):
- T4b reference, n = 4000:
  - W = n/4 = 1000 (old rule): half-window drift Cd 0.72 %, Cl 5.54 %.
    Cl fails the old tolerance max(0.5 %, 0.002).
  - W = n/2 = 2000: drift Cd 0.00350 (0.88 %), Cl 0.00073 (1.11 %).
- T4a reference, n = 3000, W = 1500: drift Cd 0.03 %, Cl 0.09 %.

New rules, oscillatory cases only (T4a, T4b, T5). T0-T3 are unchanged.
Constants are in bench/run_bench.py (STAT_*, OSC_TOL, FIELD_TOL).
- Window: W = max(1000, n/2), capped at n (n = iterations run). Was
  max(500, n/4) <= n/2. iters_to_stationary is unchanged: the run's
  window W, derived from the total n, slides over the history in steps
  of 50.
- Stationary: the two half-window means differ by at most
  max(1 % |mean|, 0.005), for Cd AND Cl. Was max(0.5 % |mean|, 0.002).
- Comparison tolerances are unchanged: Cd within max(2 %, 0.002),
  Cl within max(2 %, 0.01) of the reference window mean.
- The criterion name in the records is now
  "stationaryMean (D-042 addendum)".

Re-evaluation of the references with the new rule (read-only):

| reference | n | W | drift Cd | drift Cl | tol | stationary | iters_to_stationary | Cd mean +- std | Cl mean +- std |
|---|---|---|---|---|---|---|---|---|---|
| ref_T4a_np10 | 3000 | 1500 | 0.00011 (0.03 %) | 0.00007 (0.09 %) | 0.005 | yes | 1550 | 0.39640 +- 0.00153 | 0.07675 +- 0.00172 |
| ref_T4b_np10 | 4000 | 2000 | 0.00350 (0.88 %) | 0.00073 (1.11 %) | 0.005 | yes | 2050 | 0.39962 +- 0.00354 | 0.06568 +- 0.00322 |

Mean-field delta comparison (new):
- Averaging. The fieldAverage function object is part of the T4 and T5
  template controlDicts:
  - fields: U (mean, prime2Mean) and p (mean only), base iteration;
  - entry path functions.fieldAverage;
  - simpleFoam and coupledFoam run it from the same controlDict.
  The harness sets its timeStart per run to n - W + 1 (n = budget;
  run_bench.field_average_start), always with
  `foamDictionary -disableFunctionEntries`. The timeControl of a
  function object is active from time >= timeStart - 0.5 deltaT, so
  n - W + 1 (not n - W) averages exactly the W iterations of the force
  window. Proxy check: totalIter 150 for n = 300, timeStart 151.
- Utility applications/utilities/coupledFieldCompare:
  - Independent of libcoupledFoam. It installs to $FOAM_USER_APPBIN but
    is not in ./Allwmake yet: build it with
    `wmake applications/utilities/coupledFieldCompare`.
  - Arguments: -reference <case> -time <t> -referenceTime <t>
    -Uinf <U> -pref <p> [-parallel].
  - In parallel both cases must be decomposed identically; every rank
    reads processorN of the reference.
  - Fatal if the cells per rank, the patches or the cell centres differ
    (1e-9 relative). The proxy showed why this check matters: two scotch
    decompositions of the same T2 mesh differed (10283 vs 10270 cells on
    rank 0). T4/T5 are safe because every case copies the decomposed
    cached mesh.
  - Written at -time: UMeanDelta, pMeanDelta, magUMeanDeltaRel
    (|dU|/U_inf), CpMeanDelta (dp/p_ref) and, if both cases have it,
    UPrime2MeanDelta. Boundary values are included, so the delta Cp on
    the body is visible in ParaView.
  - fieldCompare.json (all reductions global):
    - volume-weighted RMS and max of |dUMean|/U_inf and |dpMean|/p_ref;
    - the fraction of cells (and of the volume) with
      |dUMean| > 0.05 U_inf;
    - per wall patch and over all walls: the area-weighted RMS, mean and
      max |.| of dpMean/p_ref;
    - the volume RMS of 0.5 tr(dUPrime2Mean)/U_inf^2.
  - The harness takes U_inf from forceCoeffs magUInf and uses
    p_ref = 0.5 U_inf^2.
- Proposed pass criteria (the user may change them): volume RMS
  |dUMean|/U_inf <= 0.02 AND volume RMS |dpMean|/p_ref <= 0.02
  (run_bench.FIELD_TOL). tests/test_T4_motorBike.py (T5 inherits) runs
  the utility after both runs exist (mean_field_comparison) and adds the
  checks fieldU / fieldP and the JSON to the record. The validation
  table has a new column "field RMS dU / dp".
- Reference continuation: the cached T4a/T4b references were run without
  averaging and are not recomputed from scratch.
  `continue_reference_with_average(case, n_extra)` in
  tests/test_T4_motorBike.py does the following:
  1. Continues the cached simpleFoam reference from its final time t0
     for n_extra iterations through the Allrun -restart path
     (log.simpleFoam.restart).
  2. Sets endTime t0 + n_extra and writeInterval n_extra.
  3. Adds the fieldAverage entry if it is missing, with timeStart t0 + 1
     and restartOnRestart true.
  4. Writes the rank timing to timing_continuation, recorded separately
     as continuationWallSeconds / continuationCpuHours.
  5. Evaluates the forces over the continuation only (mean, std,
     half-window drift).
  6. Records the result in reference.json under "continuation".
  Defaults for n_extra = the W of the reference budget: T4a 1500,
  T4b 2000, T5 2500. The re-evaluated reference criterion keeps using
  only the original run (force history up to t0). The test triggers the
  continuation once, when a cached reference has no mean fields; a
  completed continuation is reused. It was verified on the T2 proxy only
  (+100 iterations: fresh average of exactly 100 iterations, forces
  evaluated over 301-400). The heavy continuations are left to the
  lead.
- Proxy verification (T2 copies run/fc_proxy_cf, run/fc_proxy_sf,
  300 iterations each, averaging from 151). Serial and 2-rank metrics
  agree to 7e-15 relative:
  - volume RMS |dU|/U_inf = 0.02093;
  - volume RMS |dp|/p_ref = 0.01385;
  - max |dU|/U_inf = 0.339;
  - 13.5 % of the cells above 0.05 U_inf.
  The delta fields read back with postProcess fieldMinMax, serial and
  parallel.

## D-046 - Relaxed T1/T3 criteria and a well-defined T2 ratio (user, 2026-09-22)

The user approved a slight relaxation ("die test kriterien koennen
leicht gelockert werden. sie scheinen mir zu gering"). Evidence comes
from the acceptance runs with the D-043/D-044 settings (FABLE_REVIEW 4b).

T1 - R target 1e-6 -> 1e-5 within 400 iterations:
- Measured final R: 2.7e-6 (np1) and 4.3e-6 (np4).
- dp is within 0.057 % of simpleFoam.
- Unchanged: dp tolerance 1 %, CFL >= 100 reached, no late cuts, no
  rollbacks, and the cross-rank dp check.

T3 - Cd/Cl tolerance 0.5 % -> 2 %:
- This matches the user-approved 2 % of the wake cases (D-042).
- nDynFinal == 0 is unchanged.

T2 - iteration-ratio definition fixed, no threshold changed:
- The simpleFoam reference never reaches R < 1e-5 on all of p, Ux, Uy,
  k and omega; coupledFoam gets there at iteration 471.
- The ratio was then undefined, and the test FAILED although coupledFoam
  was the one that converged.
- Now the reference's iteration count serves as a lower bound (recorded
  as iterationRatioIsLowerBound).
- Required ratio >= 2 and the xr tolerance of 2 % are unchanged.

T0 is unchanged.

### D-040 addendum - the explicit case value overrode the fix (2026-09-22)

All six case templates carried an explicit `coupled.bounds.nutMaxFactor
1e5`, which overrides the coupledDefaults value, so the D-040 fix never
reached any case run (the T3 root-cause experiment set 1e8 explicitly).
Templates now say 1e8. Found together with a harness defect:
tests/cflib/case.py set_entry rewrote fvSolution without
-disableFunctionEntries and dropped `#sinclude "relaxation"` (k/omega
unrelaxed) in every pytest case whose fvSolution it edited (T4/T5 budget
sets); now fixed for fvSolution. T4 coupledFoam results before this
commit are invalid.

### D-042 addendum 2 - coupledFoam budget 800 on the wake cases (user, 2026-09-22)

User: 2000 coupled iterations are far too many; 400-800 were the
expectation. Evidence (read-only on run/T4a_np10, window W, drift rule
of addendum 1): coupledFoam is stationary from iteration ~400-600 (drift
<= 0.1 %), window means at 600/800/1000/1500 all within Cd
0.4029-0.4035 - much calmer than simpleFoam, whose slow wake swing needs
1000+ iteration windows. coupledFoam takes far larger pseudo-time steps
(CFL 500), so its window in iterations can be much shorter.
Change: coupledFoam budget 800 on T4a/T4b/T5 (was 1500/2000/2500);
window minimum 1000 -> 300, i.e. W = max(300, n/2): coupledFoam W = 400
(fieldAverage from iteration 401); the simpleFoam reference windows are
unchanged (n/2 = 1500/2000/2500 >= 1000 anyway). T4a result so far
(1500-iteration budget) remains valid evidence; it is rerun with 800.

## D-047 - Static remediation thresholds relaxed; test limit 1.5 % (user, 2026-09-22)

On the snappyHexMesh motorBike mesh (T4a) the static remediation set
held 1.37 % of the cells, which failed the 1 % test limit. The static
set marks cells by mesh quality at start-up and demotes them to
first-order convection with a reduced local CFL.

The user decided to raise the non-orthogonality threshold to 85 degrees
(explicitly 85, not 75) and to give the other quality criteria slightly
larger tolerances, so that fewer cells enter the set.

| threshold | old | new |
|---|---|---|
| nonOrthThreshold | 70 | 85 |
| skewThreshold | 4 | 6 |
| volRatioThreshold | 20 | 30 |
| aspectThreshold | 1000 | 2000 |
| test limit MAX_STATIC_FRACTION (T4/T5) | 1 % | 1.5 % |

The new values apply in coupledDefaults.H and in every case template.
The static non-orthogonal limiter (0.2) and the static CFL factor (0.5)
are unchanged.

All T4 coupledFoam results obtained before this entry used the old
thresholds and are rerun.
## D-045 - Deep diagnostics logging, levels 0-3 (TASK 5, 2026-09-22)

User request: an optional mega-verbose log of every sub-step, so that
parameters can later be tuned to the state of the run.

**Control.** Sub-dictionary `coupled.diagnostics`:

| keyword | values | default |
|---|---|---|
| `level` | 0-3 | 0 (off) |
| `echo` | yes/no | no |
| `maxBytes` | bytes per rank | 2 GiB |
| `upLeg` | first/all | first |

Constants live in coupledDefaults.H: `diag*`, and `forceCoeffsDrift*` for
item b below.

**Output.**
- File: `<case>/diagnostics/diag.rank<N>.jsonl`, one per rank, JSON Lines.
- Line 1 is a `header` record. It holds the level, rank, restart flag and
  the list `localKeys` of rank-local fields.
- Then one `iter` record per outer iteration. Doubles are written with
  `%.17g` (max_digits10); non-finite values become `null`.
- The file is flushed once per iteration.
- If `maxBytes` is exceeded, one `truncated` record is written and logging
  stops. The run continues.
- A restart appends to the existing file.
- `echo yes` also prints the level-1 part pretty-printed into the log
  (`CFdiag|`).
- The `CF|` lines themselves are unchanged. Item b only appends fields.

**Design.**
- Class `src/io/diagnostics.{H,C}` contains a streaming JSON builder
  (`diagJson`), the writer, and the phase classifier (`diagPhase`).
- The solver hands a pointer to the Krylov solver (`blockSolver::
  setDiagnostics`). The Krylov solver passes it on to blockGAMG through
  blockGAMGPrecon. The assembler gets `setTiming(true)`.
- Level 0 guarantee (every hook is guarded at its call site with
  `active(n)`):
  - the solver and the assembler get no pointer and no timing;
  - no strings are built, no fields are gathered, no collectives are
    added, and there are no clock reads;
  - no `diagnostics/` directory is created.
- Level 0 and levels 1-3 give identical solver trajectories: the serial
  T0 `CF|` lines are equal to those of a level-0 run (tests/test_diagnostics.py).
- The GAMG per-sweep logging at level 3 splits `smooth(x, b, n)` into n
  single sweeps. This is equivalent for both smoothers.

**Phase** (TASK 5.2, field `phase`; rules checked in this order):
1. `startup` if the global blending factor beta < 1.
2. `stalled` if the minimum of R over the last 50 iterations is not below
   the minimum of all earlier iterations. In other words: no new minimum
   for 50 iterations.
3. `asymptotic` if R <= 100 residualTol, or if the Cd/Cl window spread is
   at most 10 times its tolerance.
4. `ramp` otherwise.

Nothing tunes itself from the phase yet.

**Level 1 record** (key paths):
- Top level: `iter`, `phase`, `wallTime`.
- `residuals`: R, Rraw, rU, rp, R1, normFactor, linInitial, linFinal,
  linIts, linRestarts, linConverged, linBreakdown, rho, massErrMax*,
  massErrSum*. The mass error is sum_f phi_f per cell after the flux
  update.
- `controls`:
  - CFL, CFLstart, growth, strategy, hold, nLocLim, dt{min,median,max}*;
  - eta, etaRaw, etaClip (none | fixed | noHistory | safeguard1 | etaMin
    | etaMax);
  - omega, cuts, skipStep;
  - trials[] (CFL, eta, linIts, linConverged, linFinal, failed, omega);
  - sentinel{checks, rolledBack, nRollbacks, consecutive};
  - remediation{nStat, nDyn, version};
  - anderson{enabled, status, m, nHistory, gammaNorm, maxAbsGamma,
    flushes, flushReason};
  - beta, startupDone.
- `turbulence`: k/omega {init, final, its} from the native solver
  performance, nBoundK*, nBoundOmega*, nNutCapped, nClamped. The two bound
  counts are the cells below kMin / omegaMin before `applyBounds`. Bounding
  inside the turbulence model is not counted.
- `timings`*:
  - tAsm, split into tMomentumOps, tBoundary, tContinuity and tRhieChow;
  - tFlux;
  - tSolve, split into tPrecSetup, tPrecApply and tKrylov (= the rest,
    including Amul), with nPrecSetup and nPrecApply;
  - tTurb, tIter, tWall (the CF| value);
  - tDiag: all diagnostics work of the iteration, plus writing the
    previous record (tDiagWritePrev).
- `gamg`: nLevels, Cop, cycle, nPostSweeps, setups, reagglomerated,
  hierarchyVersion. `hierarchy` holds cellsPerLevel, ranksPerLevel, ratios,
  Cop, mergeLevels and denseCoarsest when the hierarchy changed, otherwise
  "unchanged".
- `memory`*: rssKB, peakRssKB.
- `forces`: see item b.

Keys marked * are rank-local; `bench/diag_tools.load` combines them
across ranks.

**Level 2** (key `linear`, one object per solve attempt):
- solver, initial, final, its, restarts, converged, rho;
- rhoOpt of the first preconditioner application, taken from the
  FGMRES/GMRES values of that pass at no extra cost;
- krylov[]: the residual estimate of every Krylov iteration;
- trueResidualAtRestart[]: the double-precision true residual at every
  restart / end of cycle;
- precon[] per application:
  - levels[] {l, pre[||b||, ||r||], post[before, after]};
  - coarse[] {its, res} or {dense};
  - k[] {l, r0, r1, threshold, second, a1, a2}.

The down-leg norms of every application are computed inside the cycle's
own residual: one pass that is bit-identical to the plain residual, plus
one reduction. The up-leg norms cost two extra residuals per level visit.
With `upLeg first` (default) they are only logged for the first
application of each solve (`post: null` otherwise). With every
application, the level-2 cost on T0 was 13-18 % of the iteration time
(self-measured tDiag), too close to the 15 % gate. `upLeg all`, and every
level-3 run, log every application.

**Level 3:**
- `gamgSetup`* {l, rows, dominanceMin, dominanceMedian,
  rowsWithoutOffDiagonal} after every GAMG update. Dominance =
  |a_ii| / sum_{j != i} |a_ij| per scalar row, over the diagonal and
  internal-face blocks.
- preSweeps/postSweeps norm arrays in levels[].
- trials[] gain violU/violP* (cells over the fU/fp limit at the full step)
  and violUAtOmega/violPAtOmega*.
- `bcFlips`* {patch: nFaces, nOutflow, nFlips}: faces whose sign of phi
  changed since the previous iteration, on patches whose U or p boundary
  condition is a mixedFvPatchField (inletOutlet, freestream*). This is the
  T3 chatter probe.
- `andersonInternals` {conditionEstimate = max|R_jj| / min|R_jj|,
  absGamma[], maxAlpha, maxAlphaClip, nSkippedTotal}.

**Extra parallel reductions:**
- Level 1: none. Global values are the ones the solver has already
  reduced; rank-local values stay rank-local.
- Level 2: one 2-value sum per level visit (down leg). On the up-leg
  applications, add two single-value sums per visit. Both use the level's
  communicator.
- Level 3: one sum per smoother sweep. Operator stats, violations and
  flips are rank-local.

**b. Force-coefficient window statistics** (lead/user scope extension:
"stop some cases based on the RMS of the loads, or at least monitor and
display the RMS"):
- convergenceMonitor keeps, for Cd, Cl and Cm (= CmPitch, which forceCoeffs
  always provides), over the last `forceCoeffsRmsWindow` samples (default:
  forceCoeffsWindow):
  - the mean;
  - the RMS fluctuation (population standard deviation);
  - the half-window mean drift |m_old - m_new|.
- The sums are shifted by the first window sample, so a constant sequence
  gives RMS and drift of exactly 0.
- `CF|` gets `CdMean= CdRms= ClMean= ClRms= CmMean= CmRms=` appended after
  `Cl=`. Existing fields and their order are unchanged.
- The diagnostics record gets `forces`. summary.json gets `forceStats`
  (final values) and `forceHistory` (per-iteration arrays).
- Optional stop rule, off by default (`forceCoeffsDriftTol` 0): drift <=
  max(forceCoeffsDriftTol |mean|, forceCoeffsDriftAbs = 0.005) for every
  coefficient over a full window. This is the D-042 rule. It is an
  alternative way to satisfy criterion (ii) and combines with mode any|all
  as in D-005.
- Restart: the last max(window, rmsWindow) samples of Cd, Cl and Cm are in
  coupledState. The T3 split run continues the window: n = 30 after the
  restart, from 29 restored samples plus 1 new one.
- Found on the way, fixed in `stateDict()` and convergenceMonitor: the
  restart scalars (CFL, R1, Uref, pref, PTC/EW state, Cd/Cl window) were
  rounded to 6 digits. Dictionary entries are formatted with
  IOstream::defaultPrecision() when they are created, before coupledState
  writes with max_digits10 (D-034). The state is now built with
  max_digits10.
- Test-convergenceMonitor checks, on synthetic sequences:
  - a constant: RMS and drift exactly 0;
  - a sinusoid: RMS = A/sqrt(2); the drift rule converges only when
    enabled;
  - a ramp: drift = slope W/2, no convergence;
  - the restart round trip is bit-identical.

**Measured** (T0 Re100 np1 = 57 iterations, T1 np1 = 100 iterations):

| run | level 1 | level 2 | level 3 |
|---|---|---|---|
| T0 file [B] | 109 083 | 637 610 | 1 049 142 |
| T1 file [B] | 201 612 | 349 215 | 541 380 |
| T0 bytes/iteration | 1.9 k | 11.2 k | 18.4 k |
| T0 tDiag / (tIter - tDiag) | 0.3 % | 2.4 % | ~25 % |

The last row is provisional evidence: the solver's own measurement of the
diagnostics cost, from runs pinned to a free core. It is robust against
the load of other jobs.

A wall-clock check on a lightly loaded machine (load ~2, before the
level-2 fusion and `upLeg first`) gave these medians of 3:
- pre-TASK-5 build: 6.256 s
- level 0: 6.273 s (+0.3 %)
- level 1: 6.359 s (+1.4 % vs level 0)
- level 2: 8.338 s (+33 %). This measurement led to `upLeg first`.

Measurements taken during the 10-rank T4a/T4b runs varied by +-30 %
between identical runs, so they cannot be used for the gates.

**PENDING (next quiet window, lead):**
- the wall-clock medians of the final build (gates: level 0 +-5 %,
  level 1 <= 2 %, level 2 <= 15 %);
- test_diagnostics np4;
- the np4 unit battery.

Run them unchanged with
`/home/jonas/bin/cfenv sys bash bench/task5_acceptance.sh`. It writes
run/task5_acceptance.log and run/overhead.txt.

Already passed on the final build, serially:
- test_diagnostics np1: level-3 records, level 0 writes nothing, the CF|
  lines are identical to level 0, T3 force statistics including a
  restart;
- Test-convergenceMonitor;
- test_env (Test-precision, all apps);
- Test-block4Ops;
- Test-blockFGMRES.

An earlier build of this branch also passed test_diagnostics np4 (quiet
window, before item b).

**Known issue found on the way (not caused by TASK 5):** T0 np4 is not
reproducible from run to run. In three identical runs of the pre-TASK-5
main-install binary, the iteration-1 linIters were 7/9/9 and rho
differed. The level-0 equality check of test_diagnostics is therefore
asserted serially only.

## D-051 - Amendment C design exclusions (C8, 2026-09-22)

The following were considered in an external coupled-solver comparison
and are deliberately not adopted:
- **No solve-on-variables form.** The increment form A dx = b - A x is
  required for the line search (7.2), the Eisenstat-Walker adaptive
  tolerance (B2) and the rollback (9.3) (sections 5.7, 7.2, 9.3).
- **No multiplicative relaxation stack.** The local pseudo-time step
  replaces explicit under-relaxation of the coupled system. The (2 - r)
  approximation is numerically wrong for r < 0.8 and is not adopted.
- **Static-only remediation is insufficient.** The dynamic set (8.2)
  stays.
- **No post-solve pressure shift for closed domains.** Row replacement
  (D-021) stays, because it avoids a singular system.

Also recorded: amendment C1 proposes lowering nonOrthThreshold from 70 to
65. The user set 85 on the same day (D-047), and that later, explicit
user decision stands. 65 and 60 are benchmark variants (C7) only, until
the user decides otherwise.

### D-042 addendum 3 - T5 simpleFoam reference budget 2000 (user, 2026-09-22)

The user considers 5000 simpleFoam iterations unrealistic for practice
("nobody runs that many"). The T5 Ahmed reference (coarse mesh,
1.04 M cells) therefore gets a budget of 2000 for now, with window
W = 1000. The budget is extended on the user's request if the reference
does not become stationary. coupledFoam stays at 800. The T4a and T4b
references keep their completed runs (3000 and 4000 iterations).

## D-048 - Start-up control: hybrid beta ramp (user decision, 2026-09-22)

Default startupMode hybrid:
- Upwind (beta 0) first, then a linear ramp of beta from 0 to 1 over 20
  iterations (startupRampLength). The ramp starts at iteration 10
  (startupRampStart) at the latest, or earlier once R/R1 < startupSwitchR
  (1e-2). Beta therefore reaches 1 by iteration 30 at the latest.
- The delay while R still drops fast (startupRampStartMax > 10) and the
  stagnation trigger (startupStagnationTrigger) are OFF by default: both
  measured worse.
- A developed or mapped start is detected when max(rU, rp) of
  iteration 1 is below startupDevelopedTol 1e-2 with a non-uniform field;
  beta is then 1 from iteration 1. This avoids first-order upwind on
  mapped solutions, where the relative residual never drops (user
  concern).
- startupMode upwind reproduces the old behaviour bit-identically;
  startupMode none sets beta 1 from iteration 1. Restart continues the
  ramp state.

Evidence on T1 serial (iterations to R < 1e-5):

| variant | iterations | wall |
|---|---|---|
| old upwind jump at 50 | 290 | 26.0 s |
| plain ramp 10 -> 30 | 348 | 36.2 s |
| delayed ramp | 375 | 41.1 s |
| stagnation-triggered ramp | fired at iteration 7, missed R < 1e-5 in 400 | - |

T0 Re1000 np4 did not converge within 300 iterations with the ramp in
one run; that count is chaotic (FABLE item 1). The final choice between
the ramp and the upwind jump follows the T4a comparison (pending); the
user keeps the ramp as the default meanwhile.

## D-049 - Block ILU0 pivot-growth guard (T4b collapse root cause)

Root cause of the T4b linear-solve collapse at outer iteration 65, proven
by A/B on the exact failing system:
- On a few pressure/continuity rows, mostly on GAMG level 1, the
  ILU-modified 4x4 diagonal block loses its weight. Its inverse grows up
  to 223x the inverse of the original block.
- The damped smoother then amplifies the residual by 1e3-1e5 per call,
  FGMRES stagnates, and CFL cuts cannot help because the pressure
  couplings cause it.

Guard: when the modified inverse exceeds pivotGrowthLimit (default 20)
times the inverse of the original block, the cell falls back to the
original block. The number of fallbacks is logged as nPivFb (CF| line)
and pivotFallbacks (summary JSON).

Why 20 and not 5: with 5, T2 replaced 7-10 % of its pivots (1.94 M
fallbacks), needed 25 instead of 12 linear iterations per solve and 60
CFL cuts instead of 1, and ran 3x longer. From 20 upward T2 is bitwise
unchanged, and the degenerate T4b pivots (growth 31-223) are still
caught.

A/B on the failing T4b system: 7 FGMRES iterations with the guard
instead of 3 x 400 failed. T4a with the guard: PASS.

## D-050 - Uref from the non-coupled boundary values

Uref (lineSearch::setReference) was max |U| over the field, including
the potential-flow peak: 176.5 m/s on T4b and 85.9 on T4a, against
20 m/s inflow. That effectively disabled the D-044 dynamic set on T4
(nDyn was always 0) and loosened the line search, the local limit and
the sentinel.

Now: coupled.Uref boundary (default) = the maximum over the non-coupled
boundary values; field = the old behaviour; or an explicit value.
Result on T4a: Uref 20.17 m/s, and nDyn > 0 (2523 cells at iteration 1,
then 250-330). T1 and T2 are bitwise unchanged.

Also fixed: the abort fields now land in `<iter>_lastValid` through
explicit writes. Before, v2606 redirected the write into the current
time directory. Checked serial and np4.

## D-053 - Tensorial Rhie-Chow diffusivity (amendment C2)

coupled.rhieChow.tensorial (default yes):
- D_P = V_P A_P^-1, computed in double from the 3x3 momentum diagonal
  block after PTC, including MRF Coriolis.
- Singular blocks (|det A| < 1e-12 |tr A/3|^3) fall back to a Jacobi-SVD
  pseudo-inverse (singular values <= 1e-6 sigma_max are dropped). Their
  count nPseudoInverse fails T0-T3 if it is above 0.
- D_f = n . interp(D) . n enters g_f and q_f unchanged (D-035). The scalar
  D = V/abar remains the row-scaling reference.
- tensorial no reproduces the previous scalar path exactly.

Measured on T0-T4a:
- D is isotropic except in cells at slip/symmetry patches: 0-0.14 % of
  the cells, anisotropy <= 1.4 %.
- Accuracy is unchanged, and the cost is 1-2 % of wall time.
- It is expected to matter on MRF cases (rotating F1 wheels) and at
  slip/symmetry patches.
- Kept as the spec default; revisit after the E-rcScalar benchmark on
  T4b/T5.

C4: `agglomerator algebraicPair` is the amendment-C name for
agglomerationWeights pressure (the p-p block magnitudes). Giving it
together with a different explicit agglomerationWeights is an error.
The default follows the E-algPair benchmark on T5.

## D-054 - Amendment C1/C5/C6 (static topological criteria, built-in zones, diagnostic fields)

C1 - static topological criteria (src/control/staticCriteria), OR-ed
into the static set:
- wallStarved: 3D only; a wall cell with at most 2 internal faces
  (internal, processor and cyclic faces count; cyclicAMI faces do not, so
  the set does not depend on the decomposition).
- procAMI: a cell with both a processor face and a cyclicAMI face (this
  one depends on the decomposition by construction).
- volumeJump: volume ratio > 50 across an internal, processor or cyclic
  face.

remediationFlag bits: 1 static, 2 dynamic, 4 wallStarved, 8 procAMI,
16 volumeJump. All criteria are on by default. nonOrthThreshold stays
85 (D-047/D-051).

Measured: T0-T3 add 0 cells and give identical results. T4a adds 12
wallStarved cells (static set 2655 -> 2667, 0.75 %), the same on 1 and
4 ranks. With volRatioThreshold 30, volumeJump adds no cells but keeps
its own bit and count.

C5: built-in names in zonal.zones (_wallCells, _procCells, _amiCells,
_procAMICells, _wallStarved, _remediationStatic; literal names only).
Factors of overlapping entries multiply (D-027). Regex patch sets
already worked.

C6: localDt, localCFL, cflFactorEff and betaEff (those of the last
assembled trial step) are written at every write time and into
<iter>_lastValid on abort, as extrapolatedCalculated fields
(diagnosticFields.enabled, default yes). USFD follows once SFD (C3) is
merged.

## D-052 - Selective Frequency Damping (amendment C3; proposed, osc-memory, 2026-09-22)

Implemented as specified in C3 (`coupled.sfd`, src/control/sfdControl.{H,C}):
chi* = chi Uref/Lref, Delta* = Delta Lref/Uref; filtered velocity USFD
(implicit-Euler low pass with the local dt_P of the accepted assembly, after
the local limit); momentum rows a_P += chi* V_P, b_P += chi* V_P Ubar_P
before PTC; Ubar reset to U with every Anderson flush (incl. rollback);
SFD-off after nHold (= ptc.nHold) accepted iterations with
R < deactivateBelowR; restart state sfdOn/sfdQuiet/sfdInitialised + USFD;
C6 hook diagnosticFields::setSFD while active; level-1 diagnostics
controls.sfd {active, chiStar, resets, maxDev (rank-local)}.

Choices the spec leaves open:
- The SFD term is excluded from the Rhie-Chow D (scalar abar and the
  tensorial A_P of D-053). Otherwise a state converged with SFD active
  carries a chi-dependent pressure dissipation: with the tensorial path it
  moved the converged T3-SST Cl by 2 % (-0.7 % vs -2.8 %).
- Activation after the start-up phase (`afterStartup yes`, default).
  Activated in iteration 1 the forcing towards the potentialFoam field cut
  the T3 CFL to 1 and held R ~ 4 for ~500 iterations; T2 aborted via B4
  (pre-merge build without the ILU0 guard).
- `resetOnFlush` (default yes = spec) and `startIter` (earliest activation,
  default 0) as variants for studies.

Evidence (serial, merged build = main 05f7467 + osc, 3000 its max):
- T3-SST: the only mechanism that removes the limit cycle. Converges in
  628 its / 65 s (R 9.6e-7, nDyn 0): Cd 0.0910 (+0.4 %), Cl 0.2460
  (-2.8 %) vs simpleFoam 0.0907 / 0.253. Cl amplitude over the last 200
  iterations 0.005 (without SFD: 1.04, Cd 0.123 +- 0.114). At convergence
  max|U - Ubar| = 1.1e-6 Uref: the forcing has vanished (unbiased).
  "enabled no" does not converge, so the C3 1e-3 comparison is made
  against the reference instead.
- The spec deactivation (1e-4) is harmful on T3: the normalised R is
  below 1e-4 already during the start-up ramp (SFD-off at iteration 39)
  and the limit cycle returns (identical to no SFD). Pre-merge, with a later
  SFD-off (iteration 775), R jumped from 9e-5 to 0.18 within 100 its.
  The case templates therefore use deactivateBelowR 0 (never off); the
  coupledDefaults value stays at the spec's 1e-4 (user decision).
- T3-GEKO: T3 is multi-stable in coupledFoam. SFD from the end of the
  start-up converges (487 its) to a stalled branch (Cd 0.079, Cl 0.21 vs
  reference 0.0343 / 0.838). With startIter 300, or with the implicit local
  limit (D-055), it converges to Cd 0.0363 (+6 %), Cl 0.788 (-6 %) in
  482 its.
- T0: not worse (Re100 64 vs 65 its, Re1000 79 vs 99). T1: first R < 1e-5
  at 396 vs 387 (limit 400). T2: first R < 1e-5 at 454 vs 497 but late
  bursts, final R 1.7e-4 vs 8.5e-6.

Decision: default `enabled no` (T3-GEKO lands on the stalled branch, and
the T2 end state is worse). The T3 template enables it together with the
implicit local limit (D-055): both models converge with nDyn 0.
Recommended keywords when enabling: deactivateBelowR 0, Lref 1. With
Lref = chord (35 m, chi* 0.37 1/s) the damping was too weak on T3.

## D-055 - T3 bursts: root cause, local-limit and dynamic-set memory (proposed, osc-memory, 2026-09-22)

Root cause of the T3-GEKO bursts:
- The explicit solution-limited local CFL of 7.3, dU_P = |r_P| dt_P/V_P,
  overestimates the update of stiff cells by 1 + a_P dt_P/V_P. At CFL 500
  it limits 20-2400 cells even in the nearly converged state. That makes
  the pseudo-time diagonal a non-smooth function of the residual
  (V_P/dt_P ~ |r_P|).
- A burst starts as a mode growing about 1.25x per iteration from
  R ~ 2e-4, at global CFL 500, omega 1, no cuts. nLocLim grows with it
  (23 -> 1400), and the line search stops it at omega ~0.8.
- Evidence (3000 its, iterations with R > 3e-3 after iteration 200):

| variant | iterations R > 3e-3 |
|---|---|
| explicit limiter (baseline) | 761 |
| localLimit off | 20 |
| implicit estimate (limiter practically never fires) | 9 |
| memory with permanent throttling | converges at 1193 |
| fLoc 0.2 (more limiting) | 2300 |

  - The implicit-estimate row counts 8 on the merged build.
- Second mechanism: the dynamic set releases 3-7 leading-edge cells
  (marked by cSpike) after nQuietIters. Each release switches beta 0 -> 1
  and the CFL factor 0.1 -> 1 at once, which gives a sawtooth with period
  about 24 and an R floor of about 2e-4.
- Not involved:
  - ILU0 pivot growth: the D-049 guard never fired.
  - Uref: 25.905 is the boundary maximum as well, so the T3 trajectory is
    bit-identical with and without the Uref change.
  - Freestream flips: 0.
  - k bounding: fires in every iteration, also in converged states.
  - Eisenstat-Walker: eta = 0.5 throughout.
  - Global CFL, line-search cuts and sentinel: CFL stays at 500, there
    are 0 cuts and 0 rollbacks during bursts. Their release after a cut
    is already gradual: a hold of nHold, then growth of at most betaMax
    per iteration.
- The T3-SST "limit cycle" is the iteration alternating between an
  attached branch (Cd ~0.032, Cl ~0.9-1.0) and a stalled branch
  (Cd ~0.09, Cl ~0.25; the simpleFoam SST reference sits on this one).
  The Cl amplitude over 200 iterations is 1.0.

Implemented (all default off, "off" = bit-identical old behaviour):
- `localLimit.implicit` (no): Jacobi estimate
  dU_P = |r_P|/(V_P/dt_P + a_P), with a_P the mean momentum diagonal
  without PTC (coupledAssembler::momentumDiag).
- `localLimit.memory` (no), `localRecovery` 1.5, `localHold` 0,
  `localStickyAfter` 0: a per-cell dt factor f_P in (0, 1].
  - When the cell is limited, f_P is cut to what the check needs. Every
    line-search trial starts from the factor at the start of the
    iteration; at most one limit event is counted per iteration.
  - f_P recovers by x1.5 per iteration after the hold.
  - After localStickyAfter events the factor is kept permanently.
- `remediation.dynamic.stickyAfter` (0): a cell that enters the set for
  the N-th time stays in it.
- `remediation.dynamic.releaseIters` (0): release ramp. Over N iterations
  beta goes linearly 0 -> 1 and the CFL factor geometrically
  cflFactor -> 1.
- Restart state: sparse lists localLimit*, dynamicEntry*, dynamicRamp*.
  Diagnostics:
  - level 1: nLocThrottled, nLocSticky, nDynSticky, nDynRamping;
  - level 3: the top-20 locally limited cells and the dynamic-set members
    with their positions.

Evidence on the merged build (T0 Re100 / Re1000 iterations; T1 and T2:
first iteration with R < 1e-5; T3-GEKO: iterations with R > 3e-3 after
iteration 200):

| mechanism | T0 | T1 | T2 | T3-GEKO | result |
|---|---|---|---|---|---|
| stickyAfter 3 | bit-identical | bit-identical | bit-identical | converges at 404 its: Cd -6.4 %, Cl +6.9 % | opt-in |
| memory 1.5 | 66 / 63 | 394 | 561, final 3.7e-5 (worse) | 157 over (441) | opt-in |
| release ramp 10 | bit-identical | 378 | 466, but late bursts (final 5.4e-3) | 125 over | opt-in |
| implicit | 66 / 116 | 390 | 496 | 8 over | opt-in, T3 template yes |
| baseline | 65 / 99 | 387 | 497 | 441 over, not converged | - |

Notes on the table:
- stickyAfter 3 keeps 3 first-order leading-edge cells, so nDyn == 0 at
  convergence (spec 13) fails by construction. On T3-SST it converges to
  the wrong (attached) branch: Cd -65 %, Cl +288 %.
- Memory 1.5: T3-SST is unchanged.
- Implicit: T3-SST still limit-cycles. T4/T5 have not been tested.

Defaults: all opt-in. Fable's WIP default stickyAfter 3 is reverted to 0.
The T3 template uses implicit yes together with SFD (D-052).

## D-058 - SFD default never off; T3 template with SFD; T3 tolerance 5 % (user, 2026-09-22)

The user took three decisions on the evidence of D-052/D-055 (branch
osc-memory, merged):

1. **SFD never switches off by default.** coupled.sfd.deactivateBelowR
   defaults to 0 (coupledDefaults::sfdDeactivateBelowR). With the spec
   value of 1e-4, SFD switched off at iteration 39 of T3-SST (R drops
   below 1e-4 briefly during the start-up) and the limit cycle returned.
   Keeping SFD on is harmless: the forcing -chi*(U - Ubar) vanishes at
   convergence (|U - Ubar| = 1.1e-6 Uref measured), so the converged
   solution is unbiased. SFD itself stays OFF by default and is enabled
   per case.
2. **The T3 template keeps the implicit local-limit estimate and SFD
   enabled.** Both turbulence models then converge without dynamic-set
   cells: SST in 628 iterations, Cd +0.4 % and Cl -2.8 % vs simpleFoam;
   GEKO about +-6 %. This is how a user would set up such a case.
   T3 has several steady branches in coupledFoam (attached and stalled),
   and simpleFoam's heavy relaxation lands on a different branch, so the
   remaining offsets reflect the near-stall physics of this case (D-055).
3. **T3 Cd/Cl tolerance 2 % -> 5 %** (tests/test_T3_airFoil.py; history:
   0.5 % spec -> 2 % D-046 -> 5 % D-058). User-approved relaxation, not an
   agent loosening.

## D-059 - 12-hour budget of the final re-run (user, 2026-09-22)

The full final re-run (TASK 6) was estimated at 18-24 h, dominated by the
heavy benchmark (simpleFoam A/B on T4b and T5 alone about 4-5 h) and by
strong scaling on T4b (serial run several hours). The user asked for about
12 h and decided:

- **No SMT**: all timing runs on physical cores only (max 16 ranks).
- **One repetition** of the heavy benchmark cases (T1-T3: 3 repetitions).
- **T-scaling on T4a** instead of T4b, ranks 1, 2, 4, 8, 12, 16, 150
  iterations per run with the first 50 excluded
  (CF_SCALING_MESH=a, CF_SCALING_ITERS=150; results/tests/T_scaling_T4a.json).
  The pass criterion (coupledFoam efficiency >= 0.8 x simpleFoam) is
  unchanged.
- **Heavy benchmark scope**: T4a runs every configuration (A-H and the C7
  E variants, --no-scope). T4b and T5 run configuration E only; their
  simpleFoam comparison uses the simpleFoam references of the T4/T5 tests
  (run/ref_*, same system simpleFoam and tutorial settings, unchanged by
  coupledFoam commits), which the final re-run keeps.

## D-056 - Rollback restores k, omega, nut verbatim (G-lookup crash, 2026-09-22)

Symptom: T4b np10 (main 05f7467) aborted after iteration 5 on all ranks
with "failed lookup of kOmegaSST:G"; the same FATAL in the osc T3
experiments (fresh starts and restarts, kOmegaSST and GEKO).

Root cause: sentinel::restore() (the rollback of spec 9.3) called
omega.correctBoundaryConditions() outside turbulence->correct(). The
native omegaWallFunctionFvPatchScalarField::updateCoeffs() (and
epsilonWallFunction) looks up the production field <model>:G, which
kOmegaSST, GEKO and kEpsilon register only temporarily inside correct().
Every rollback of a wall-function case therefore crashed; the log shows
no "sentinel rollback" line because the message is printed after
restore(). Checked and NOT involved: bound() in applyBounds (v2606
bound.C assigns the boundary values and evaluates only coupled patches -
T4a/T4b bound omega every run without a crash), the restart read path
(fields are read, the model's validate() only corrects nut; the osc
restart crashes happen after the first linear solve and before the
turbulence solve, i.e. at the first sentinel check), nut wall functions
(no G lookup), diagnosticFields, lastValid writes (copies, no evaluate).

Fix: store() keeps the boundary values of k, omega and nut as well;
restore() assigns internal and boundary values with a forced assignment
(==), no updateCoeffs()/evaluate(). The stored values are those the model
itself evaluated at the end of the previous accepted iteration, so the
restored state is exactly that state (processor patches included, all
ranks restore together). U and p keep correctBoundaryConditions().

Verification: T3-SST forced rollback (sentinel.UFactor 1.2): crash at
the first rollback (iteration 11) before, controlled rollbacks after
(recoveries, then the 9.3 abort with <n>_lastValid when the limit stays
violated). On the merged build (SFD + implicit local limit): T3-SST and
T3-GEKO forced-rollback runs and restart-with-rollback runs - no G
failure. The restart round trip itself is bitwise identical between the
fixed and the unfixed build; it is not exact run-to-run on T3 (open
point, pre-existing, FABLE_REVIEW 6).

## D-057 - Potential-flow start: start-up references (2026-09-22)

Symptom: with Uref from the boundary (D-050) T4a and T4b iterated at
CFL 1 with omega = omegaMin (0.1) and 3 CFL cuts in every iteration of
the start (T4a 22-25 iterations, CFL 500 only at iteration ~480; T4b
every iteration until the D-056 crash). Evidence (T4a, diagnostics level
3 and field analysis of the potentialFoam start):
- The line search is PRESSURE-limited: rank 0 counted ~300 cells over
  fp pref against ~20 over fU Uref, and omega hardly depends on CFL
  (0.039 at CFL 5, 0.035 at CFL 1): the pressure has no pseudo-time term,
  so CFL cuts cannot reduce dp. max |dp| ~ 14 pref per step vs the
  allowance fp pref = 0.5 pref.
- The cells that demand it are the singular potential-flow peaks: 6 of
  the 7 cells with |dp| > 10 pref (per unit omega) have |U0| > 3 U_inf;
  the largest |dU| (3-7 U_inf) sit at |U0| = 2.8-4.3 U_inf. The potential
  start has 18 cells above 3 U_inf on T4a (peak 4.3), 165 on T4b (peak
  8.75, p down to -36 pref); the developed simpleFoam fields have max
  |U| 1.5 U_inf (T4a).
- A CFL cut at CFLmin repeated the identical solve (3 solves per
  iteration, 2.3-2.6 s instead of 0.7 s on T4a).

Variants (keywords, constants in coupledDefaults.H):
- (a) coupled.potentialClip f: clip |U| of the potential start to f Uref,
  p shifted by the lost dynamic head. T4a with f = 3 (18 cells) stays
  crippled (omegaMin 25 iterations, CFL 5.9 at 150); on T4b (f = 4, 33
  cells) it worsens the first steps (omega 0.25 instead of 0.55). Default
  0 (off).
- (b) coupled.startupReference ramp: Uref and Ustep blend from Ufield0
  (max |U| of the initial field) to their frozen values along the beta
  ramp. T4a 150 its: CFL 500 at 100, no cut.
- (c) coupled.startupReference exclude: while the start-up runs (beta < 1,
  i.e. at most 30 iterations, D-048) and only if Ufield0 exceeds the step
  scale: step scale Ufield0, no dynamic-set marking (rollback marking
  stays). T4a 150 its: CFL 500 at 100; T4b 60 its: no cut, CFL 321 at 60,
  R 1.24e-2 (best of all T4b screens).
- coupled.UrefStep reference | fieldCapped | field | <value> (+
  UrefStepCap, default cU = 4): the step scale of the whole run
  (line search, offending cells, local CFL limit), separate from the
  classification Uref (dynamic set, clipping, sentinel, Anderson norm,
  SFD). fieldCapped = max(Uref, min(Ufield0, 4 Uref)). T4a 800 its:
  stationary at 500 (CFL 500 at 97); T4b 60 its: omega < omegaMin at
  iteration 7 (cap binding at the 8.75 U_inf peaks) -> 2 cuts + hold ->
  CFL stalls at 6.6-9.5 (R 1.26e-2); uncapped field: no cut, CFL 237.
- A cut at CFLmin is no longer repeated (hold as before, one solve).

Decision (default): startupReference exclude, UrefStep reference,
potentialClip 0. Reasons: it removes the crippled start on T4a and T4b,
touches only the start-up of starts whose initial field exceeds the
reference (a no-op on T0-T3: Ufield0 <= Uref there), and leaves D-044 and
D-050 untouched in the developed solution (dynamic set cU 4 U_inf / cp 15
pref / cSpike 0.5 U_inf, sentinel 10 U_inf, line search fU/fp against
U_inf). fieldCapped is kept as an option (faster developed phase is not
proven: T4a 500 vs 450 stationarity; it fails the T4b start alone).

Results T4a np10, 800 iterations, D-042 criterion (W = 400; "under
load": osc serial jobs and my serial regression ran alongside, the wall
and CPU numbers are validation, not the timing comparison):

| run (build) | Uref / step / start-up | CFL 500 at | stationary at | wall / CPU-h to stationarity | total wall / CPU-h | Cd / Cl |
|---|---|---|---|---|---|---|
| R0 (main + D-056) | boundary / reference / none | 482 | 600 | 605 s / 1.68 | 789 s / 2.19 | 0.4004 / 0.0630 |
| R1 | field (all) | 75 | 450 | 521 s / 1.45 | 927 s / 2.57 | 0.4046 / 0.0645 |
| R2 | field + pivotGrowthLimit 5 | 88 | 450 | 526 s / 1.46 | 921 s / 2.56 | 0.4032 / 0.0645 |
| S (unbound, heavy load) | boundary / fieldCapped / none | 97 | 500 | (1032 s / 2.87) | (1418 s / 3.94) | 0.4005 / 0.0604 |
| pytest T4a (merged, default) | boundary / reference / exclude | 87 | 450 | 401 s / 1.11 | 712 s / 1.98 | 0.4026 / 0.0630 |
| earlier T4a_np10_fix (before amendment C) | field, pivot 5, delayed ramp | 76 | 450 | 422 s / 1.17 | 756 s / 2.10 | 0.4032 / 0.0794 |

simpleFoam reference: stationary at 1550, 783 s / 2.17 CPU-h; Cd 0.3964,
Cl 0.0768. Decomposition of the 25 % (lead): the D-050 step limits cost
150 iterations to stationarity (R0 600 vs R1 450); the pivot limit 5 vs
20 has no effect (R1/R2 450/450, 521/526 s); the remaining wall-clock gap
of R1 vs the earlier fix run at the same 450 iterations is per-iteration
cost (1.16 vs 0.95 s/it: machine load and/or amendment C, not isolated).
The default recovers it (450 iterations, 401 s).

T4b np10 (1.70 M cells, pytest, merged build, default): the run that
aborted after 5 iterations on main now runs its full 800-iteration
budget without abort and without a single rollback. Start: omega
0.73-1 from iteration 1, no CFL cut, CFL 500 from iteration ~100, R 8e-3
in the plateau (main: CFL 1, omega 0.1, 3 cuts and nDyn 8031-8351 in
every iteration, 10 s/iteration, then the D-056 crash). D-042 evaluation
(W = 400): stationary at 450, Cd 0.4078 vs 0.3996 (+2.04 %, tolerance
2 % - marginal FAIL by 0.0002 absolute), Cl 0.0655 vs 0.0657 (PASS),
mean-field RMS |dUMean|/U_inf 0.0061 and |dpMean|/p_ref 0.0011 (limit
0.02 each), static set 0.49 %, rollbacks 0. Wall and CPU under load (osc
serial jobs and the serial regression ran alongside; NOT the final
timing measurement): 3757 s / 10.44 CPU-h for 800 iterations, 2117 s /
5.88 CPU-h to stationarity against the cached simpleFoam reference
5942 s / 16.50 CPU-h - speed-up 2.81 wall and 2.81 CPU.
## D-060 - User-judged convergence point of the force cases (user, 2026-09-22)

The user is not yet sure about the automatic convergence criteria
(12.3(ii) window for T3, D-042 stationary mean for T4/T5). After the final
re-run the user inspects the load histories and may name, per run, an
earlier iteration from which they consider it converged. Decision details:

- The stop criteria of the runs stay unchanged (user). A user iteration can
  therefore only move the convergence point earlier; a value beyond the last
  iteration of a run is ignored and flagged.
- File `report/user_convergence.json` (by hand), keys T3-SST, T3-GEKO, T4a,
  T4b, T5; entries "coupledFoam"/"simpleFoam" (test run and its reference)
  or a benchmark configuration letter. `bench/make_report.py` writes
  `report/user_convergence_template.json` with the automatic iterations as
  hints.
- With a user iteration N: iterations, wall time and CPU-hours to
  convergence count up to N (run_bench.to_convergence), Cd/Cl are the means
  over N..end, speed-ups and reference deviations are recomputed.
- The pass/fail of the tests is NOT changed; it keeps the automatic
  criteria. Table `convergence_choice` lists automatic and user points side
  by side; the load plots mark the user point with a solid line.
