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
