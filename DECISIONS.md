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
