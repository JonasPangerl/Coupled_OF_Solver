<!-- Plain-text extraction of SPEC_coupledFoam.html (artifact https://claude.ai/artifact/MYAB7vRgzPqrqjtWKpERpT, fetched 2026-09-21).
     The HTML file is the authoritative copy. Section 5.5 is truncated in the source itself; see DECISIONS.md D-002. -->



coupledFoam — Implementation Specification v2

Status: Reviewed and corrected version of plan_coupled_solver.md (v1). This document is the single source of truth for a one-shot implementation by an autonomous coding agent. It is written to leave no room for interpretation. Where v1 and v2 disagree, v2 wins.

Target: Block-coupled p–U solver for steady, incompressible, turbulent external aerodynamics in OpenFOAM ESI v2606, CPU-only, single-precision linear algebra (SPDP), 16 cores / 128 GB RAM, meshes up to 45 M cells.

Deliverables: (1) the solver + library, (2) a fully automated test suite, (3) a benchmark harness against native OpenFOAM solvers, (4) an English report with comparisons and plots.

0. Rules for the implementing agent

These rules are mandatory and override any implicit judgement.

Phases and gates are sequential. Do not start phase N+1 before every acceptance criterion of phase N is met and recorded in results/gates/phase_N.json. If a gate cannot be passed, stop, write results/gates/phase_N_FAILED.md with the root cause and the exact reproduction command, and do not proceed.

Never change an acceptance threshold in this document to make a gate pass. Thresholds may only be tightened.

Never silently drop a feature. If something is infeasible, record it in DECISIONS.md with the reason and the impact.

Every numerical result (residuals, force coefficients, timings, memory) is written to JSON under results/ by the harness scripts, never only to stdout.

All parameters are runtime dictionary entries with defaults defined in Section 11. No hard-coded magic numbers in source.

Compile with warnings as errors for the project code (-Wall -Wextra -Werror in Make/options, OpenFOAM headers excluded via -isystem).

Commit after every completed subsection with the message pattern phase-<N>: <what>; keep CHANGELOG.md current.

Reproducibility: every case has an Allrun, an Allclean, and produces deterministic outputs given the same rank count.

Floating-point exceptions are bugs in development builds. Every development run sets FOAM_SIGFPE=true FOAM_SETNAN=true. A run that dies from a FPE is a failed test, not an acceptable outcome.

Read before writing: before implementing any block-matrix or AMG component, read the corresponding native OpenFOAM source (src/OpenFOAM/matrices/lduMatrix/, src/OpenFOAM/matrices/LduMatrix/, GAMG/) and follow its conventions (LDU addressing, interface handling, solveScalar).

1. Review of the v1 plan — what changed and why

#
v1 statement
v2 decision
Reason

R1
Port ICSFoam as the code base
Clean-room implementation; ICSFoam is a design reference only (read it, do not depend on it)
ICSFoam targets an older ESI release, is density-based first, and stores in scalar. A port of unknown cost is not acceptable in a one-shot. Clean-room gives full control over single-precision block storage and SPDP.

R2
Selective block-AMG (SAMG, Uroić/Jasak) as primary linear solver
Additive-correction block-AMG ("block-GAMG") reusing native GAMGAgglomeration as primary; SAMG is Phase 2
Native agglomeration (faceAreaPair, processor-aware) already exists and is proven; only the coarse-operator (block Galerkin sum) and block smoothers are new. SAMG needs a full coarsening/interpolation implementation plus parallel setup — it is the biggest risk item and belongs behind a working baseline. This is the Mangani et al. (2014) approach, which already demonstrated near-linear scaling.

R3
petsc4Foam FieldSplit as reference
Optional (Phase 2), not part of the gate path
External dependency (PETSc single-precision build) must not block delivery. Native simpleFoam is the mandatory benchmark reference.

R4
Memory budget 100–115 GB at 45 M cells
Kept, but measured at every gate; the 45 M case is not part of the one-shot — the largest one-shot case is ~8 M cells
Nothing above ~8 M cells is needed to validate correctness and measure speedup; the 45 M run is a user acceptance step after delivery.

R5
AMR as Phase 2
Kept as Phase 2, and the restart format is designed AMR-capable from day 1 (Section 10)
Unchanged.

R6
Tests listed as levels
Concrete cases, commands, acceptance numbers, and a benchmark protocol (Sections 13–15)
v1 was not executable.

R7
CFL strategy "mRDM default"
Kept, with exact formulas (Section 7)
v1 named methods without defining them.

R8
Remediation cell sets
Kept, with exact detection formulas, hysteresis, and diagnostics (Section 8)
Same.

R9
FPE section
Extended with exact guard list, double-accumulated reductions, and a rollback protocol (Section 9)
Same.

R10
Ahmed body / motorBike
Added as T5/T4 with geometry generation script and mesh targets
Requested.

2. Scope

2.1 In scope (one-shot)

Library libcoupledFoam and application coupledFoam (steady, incompressible, RANS, block-coupled p–U, segregated turbulence).

Single-precision block linear algebra (solveScalar) with double-precision reductions.

Block-GAMG preconditioned block-BiCGStab / block-GMRES(10).

Pseudo-transient continuation (PTC) with adaptive CFL (mRDM, EXP), physicality line search, solution-limited local CFL.

Two-tier remediation cell sets (static mesh-quality, dynamic solution-spike).

Implicit boundary conditions, implicit MRF, cyclic/AMI as explicit-coupled (see 5.5).

FPE safety: trapping in dev builds, guards in production builds, sentinel + rollback.

Restart (exact resume of CFL state, remediation set, phi consistency).

Test suite (unit + cases T0–T5), benchmark harness, English report.

2.2 Out of scope (Phase 2, Section 17)

AMR, SAMG, PETSc reference, GPU, compressible, transient/ALE.

3. Environment and toolchain (exact)

Item
Requirement

OS
Linux x86_64 (Ubuntu 24.04 or equivalent)

Compiler
GCC ≥ 13 (WM_COMPILER=Gcc), flags for project code: -O3 -march=native -fno-math-errno (no -ffast-math: it breaks NaN/Inf detection)

MPI
OpenMPI ≥ 4.1, WM_MPLIB=SYSTEMOPENMPI

OpenFOAM
ESI v2606, built from source with WM_PRECISION_OPTION=SPDP, WM_LABEL_SIZE=32, WM_COMPILE_OPTION=Opt. A second build with WM_COMPILE_OPTION=Debug and SPDP is required for FPE-trapping test runs. Both builds live side by side (OpenFOAM-v2606-Opt, OpenFOAM-v2606-Debug).

Python
3.11+, packages: numpy, pandas, matplotlib, pytest, pyyaml, numpy-stl (for Ahmed STL generation)

Hardware
16 physical cores, 128 GB RAM. All benchmarks: mpirun --bind-to core --map-by core -np <N>; hyper-threading disabled or ranks pinned to physical cores.

Native reference solver
simpleFoam from the same v2606-SPDP-Opt build (identical precision, identical compiler flags).

Phase 0 must verify with Test-precision (write a 5-line app) that sizeof(Foam::solveScalar)==4 and sizeof(Foam::scalar)==8.

4. Repository layout

coupledFoam/
├── README.md                     # build + run quickstart
├── CHANGELOG.md
├── DECISIONS.md                  # every deviation/decision with reason
├── SPEC_coupledFoam.md           # this document (copied in)
├── Allwmake / Allwclean
├── src/
│   ├── blockMatrix/              # lib: block LDU storage, ops, interfaces
│   │   ├── blockLduMatrix4.{H,C}
│   │   ├── blockLduInterface*.{H,C}
│   │   ├── block4Ops.H          # inline 4x4 dense ops (LU, inverse, matvec) in solveScalar
│   │   └── doubleReduce.{H,C}   # double-accumulated norms / dots with MPI reduce
│   ├── blockSolvers/
│   │   ├── blockSolver.{H,C}     # base class, runtime selection
│   │   ├── blockBiCGStab.{H,C}
│   │   ├── blockGMRES.{H,C}
│   │   ├── blockGAMG.{H,C}       # additive-correction block AMG (V-cycle)
│   │   ├── blockSmoothers/{blockGaussSeidel,blockILU0}.{H,C}
│   │   └── blockPreconditioners/{blockDiagonal,blockGAMGPrecon}.{H,C}
│   ├── assembly/
│   │   ├── coupledAssembler.{H,C}  # builds the 4x4 system from fvMesh + fields
│   │   ├── convection.{H,C}        # deferred correction
│   │   ├── rhieChow.{H,C}
│   │   ├── boundaryCoupling.{H,C}  # implicit BC rows/cols
│   │   └── MRFCoupling.{H,C}
│   ├── control/
│   │   ├── ptcControl.{H,C}        # local dt, CFL strategies
│   │   ├── lineSearch.{H,C}        # physicality omega
│   │   ├── remediation.{H,C}       # static + dynamic cell sets
│   │   └── sentinel.{H,C}          # finite checks, rollback
│   └── io/
│       └── coupledState.{H,C}      # restart state dictionary
├── applications/
│   ├── coupledFoam/                # the solver
│   └── test/                       # OpenFOAM-style unit-test apps (Test-block4Ops, Test-blockMatrix, Test-blockGAMG, Test-precision)
├── cases/                          # T0..T5 case templates (Allrun/Allclean)
├── tests/                          # pytest: unit + case tests, JSON outputs
├── bench/                          # benchmark harness + report generator
│   ├── run_bench.py
│   ├── rank_wrapper.sh             # /usr/bin/time -v per rank -> RSS
│   ├── parse_forces.py
│   └── make_report.py
├── results/                        # generated JSON, figures (git-ignored except gates/)
└── report/REPORT.md                # generated

5. Numerical specification

5.1 Equations and notation

Steady incompressible RANS, kinematic pressure p (m²/s²), velocity U, effective viscosity ν_eff = ν + ν_t.

Momentum: ∇·(U U) − ∇·(ν_eff ∇U) − ∇·(ν_eff (∇U)ᵀ) + ∇p = S_MRF

Continuity: ∇·U = 0

Turbulence (k, ω, ν_t) is frozen during the coupled solve and updated segregated after it (5.8).

Notation: cell P, neighbour N, face f with area vector S_f (outward from owner), |S_f|, face distance vector d_f = x_N − x_P, interpolation weight w_f (owner weight, from mesh.weights()), face flux φ_f (m³/s), cell volume V_P.

Unknown vector per cell: x_P = (u, v, w, p)ᵀ.

5.2 Block system layout and storage

blockLduMatrix4 stores, for nCells cells and nFaces internal faces (LDU addressing from mesh.lduAddr()):

Array
Type
Size
Layout

diag_
solveScalar
16·nCells
row-major 4×4 per cell, index cell*16 + r*4 + c

upper_
solveScalar
16·nFaces
coefficient of owner-row / neighbour-column, same 4×4 layout

lower_
solveScalar
16·nFaces
neighbour-row / owner-column

source_
solveScalar
4·nCells
RHS

interface coeffs
solveScalar
16·nPatchFaces per coupled patch
as native interfaceBouCoeffs/interfaceIntCoeffs, but 4×4

Row/column meaning: rows 0–2 = momentum x,y,z; row 3 = continuity. Columns 0–2 = u,v,w; column 3 = p.

All matrix coefficients are assembled in scalar (double) locally and stored as solveScalar (float) — the conversion happens once per outer iteration at store time. Fields U, p, phi stay scalar.

5.3 Assembly (per outer iteration)

All contributions are written into the four sub-blocks. Sign convention: A x = b, with A the full block matrix.

(a) Momentum–momentum (rows 0–2, cols 0–2):

Convection, deferred correction: implicit part first-order upwind on φ_f (from previous iteration, Rhie–Chow consistent): upper[f] += max(-φ_f,0)... precisely: owner diag a_P += max(φ_f, 0), owner/neighbour coupling upper[f] = −max(−φ_f, 0), lower[f] = −max(φ_f, 0), neighbour diag a_N += max(−φ_f, 0) (standard OpenFOAM upwind fvm::div). Explicit high-order correction: b_P −= β_P · Σ_f φ_f (U_f^{HO} − U_f^{UD}) where U_f^{HO} is evaluated with the user scheme (default linearUpwind grad(U)), and β_P ∈ [0,1] is the per-cell blending factor (1 = full high order; reduced by remediation, Section 8, and by the start-up switch, 7.4). This is implemented by evaluating fvc::div(phi, U, "div(phi,U)") minus fvc::div(phi, U, "div(phi,U)_upwind") and multiplying by β.

Diffusion: standard fvm::laplacian(ν_eff, U) coefficients (scalar, applied identically to the three momentum rows on the diagonal of the 3×3 sub-block), including non-orthogonal correction as explicit source (limited by nonOrthLimiter, default 0.5). The transpose term ∇·(ν_eff (∇U)ᵀ) is explicit (fvc::div(ν_eff · dev2(T(grad(U))))).

PTC: a_P += V_P / Δt_P on all three momentum diagonals; b_P += (V_P / Δt_P) U_P^{old} (Section 5.4).

MRF: 5.6.

(b) Momentum–pressure (rows 0–2, col 3): pressure gradient, Green–Gauss face-based: contribution of face f to owner row: + w_f S_f on the owner column-3 diag entry and + (1−w_f) S_f on the upper block (owner row, neighbour column 3); for the neighbour row the signs are reversed (S_f points out of owner). Boundary faces: 5.5.

(c) Pressure–momentum (row 3, cols 0–2): divergence of interpolated velocity: face f contributes + w_f S_fᵀ to the owner diag (row 3, cols 0–2), + (1−w_f) S_fᵀ to upper (owner row 3, neighbour cols 0–2); reversed sign for neighbour row.

(d) Pressure–pressure (row 3, col 3): Rhie–Chow Laplacian: D_f = interpolate(V/ā_P)_f where ā_P is the average of the three momentum diagonal coefficients of cell P after step (a) (including PTC term, excluding the pressure column). Coefficient −D_f |S_f|² / |d_f · n_f| with the standard fvm::laplacian(D, p) assembly (non-orthogonal correction explicit). Explicit RC source on row 3: b_P^{(3)} = − Σ_f D_f [ (∇p)_f^{avg} · S_f ] where (∇p)_f^{avg} is the linear interpolate of cell gradients (the "average gradient" term of Rhie–Chow). The implicit face gradient term is what the Laplacian block represents.

(e) Flux update after solve: φ_f = (w_f U_P + (1−w_f) U_N) · S_f − D_f [ (p_N − p_P)|S_f|/|d_f·n_f| − (∇p)_f^{avg} · S_f ], then fvc::makeAbsolute where MRF applies. This φ_f is used for the next assembly, for turbulence, and is written on output.

(f) Row scaling: before the solve, scale the continuity row by s_p = 1/max(D_ref, VSMALL) with D_ref = gAverage(V/ā_P) (dimensional balancing of the block for single precision). Unscale x after the solve. Record s_p in the log.

5.4 Pseudo-transient continuation, local time step

Per cell: Δt_P = CFL · V_P / max( ½ Σ_f |φ_f| + ν_eff,P · V_P^{1/3} , VSMALL ) (convective + diffusive limit). CFL is the global adaptive value from Section 7; the local cap of 7.3 may reduce Δt_P further. No PTC term on the continuity row.

5.5 Boundary conditions (implicit treatment)

For every boundary face f of patch cell P, the block row of P receives the following. Keep exactly this table; add to DECISIONS.md if a new BC type is needed.

Patch type (U / p)
Momentum rows
Continuity row

fixedValue U / zeroGradient p (inlet, wall, MRF wall)
convection: a_P += max(φ_f,0), b_P += max(−φ_f,0) U_b; diffusion: `a_P += ν_eff
S_f

zeroGradient U / fixedValue p (outlet)
convection: a_P += φ_f (U_f = U_P); pressure gradient: b_P −= S_f p_b; diffusion: none
divergence: + S_fᵀ on own U columns; p-p: `a_P^{(3,3)} −= D_f

slip / symmetry
normal component removed via projection: assemble as fixedValue with U_b = U_P − (U_P·n)n treated implicitly by adding `ν_eff
S_f

inletOutlet U / totalPressure p
treat per face by sign of φ_f as inlet (fixedValue) or outlet (fixedValue p with p_b = p0 − ½
U_P

cyclic / cyclicAMI
explicit-coupled in v2: neighbour values lagged from previous iteration (source), full implicit block coupling across cyclic interfaces is Phase 2
same

processor
fully implicit via block interface coefficients (as native)
same

Wall-function BCs for ν_t (nutkWallFunction, nutUSpaldingWallFunction) act through ν_eff on the wall face and need no special block handling.

5.6 MRF (implicit)

Inside MRF zones, the Coriolis term 2Ω × U is added implicitly to the 3×3 momentum diagonal block as the skew-symmetric matrix [Ω]_× scaled by 2 V_P; the centrifugal term is absorbed in p (relative frame formulation as native MRF.DDt). Fluxes on MRF-zone faces are made relative before assembly (MRF.makeRelative(phi)) and absolute after the flux update. Native IOMRFZoneList is reused for zone definitions.

5.7 Outer iteration algorithm

read fields U, p, phi, k, omega, nut; read coupledState (restart) or init
if not restart and potentialInit: run potentialFoam-equivalent init (native potentialFoam call via Allrun)
build static remediation set (Section 8.1)
for iter = 1 .. maxIter:
    turbulence->correct() is NOT called here (see below)
    compute nu_eff, grad(U), grad(p)
    compute local dt (5.4) with current CFL and local caps
    assemble block system (5.3) with beta field
    scale rows (5.3f)
    solve  A dx = b - A x   (solve for the increment dx, tolerance per Section 6)
    unscale
    omega = lineSearch(dx)                       (Section 7.2)
    if omega < omegaMin and cflCuts < maxCflCuts:
        CFL *= kappa; cflCuts++; goto assemble   (repeat the step; fields unchanged)
    U += omega dx_U ; p += omega dx_p
    update phi (5.3e)
    sentinel check (Section 9.3); on failure: rollback + CFL*=0.25 + mark cells; continue
    dynamic remediation update (Section 8.2)
    turbulence->correct()  (segregated k, omega, nut update; Section 5.8)
    residual bookkeeping (Section 12), CFL update (Section 7.1)
    write if writeTime; write coupledState
    check convergence (Section 12.3); break if converged

Solving for the increment (not the full vector) is mandatory: it keeps the single-precision solve well-conditioned relative to the field magnitude.

5.8 Turbulence coupling

Turbulence models are used unmodified via incompressible::turbulenceModel (SST, GEKO, and any user model such as SST-sf). turbulence->correct() is called once per outer iteration after the coupled update.

The k and ω matrices use native fvMatrix<scalar> with the native solvers from fvSolution (default: smoothSolver symGaussSeidel, relTol 0.1), with native relaxationFactors for equations k and ω (default 0.7). The pseudo-time term for k/ω is native fvm::ddt disabled (steady); stability comes from equation relaxation.

Bounds: after correct(), enforce k ≥ kMin, ω ≥ omegaMin from fvSolution.coupled.bounds (defaults kMin 1e-12, omegaMin 1e-6) using bound().

6. Linear solver specification

All operations on solveScalar (float). All reductions (dot products, norms) accumulate in double via doubleReduce::dot() / ::norm2() and use MPI_DOUBLE in reduce. This is non-negotiable (float sums over 10⁷–10⁸ entries lose 3–4 digits and can overflow in ‖·‖²).

6.1 Runtime selection (fvSolution.solvers.coupled)

coupled
{
    solver          blockBiCGStab;   // blockBiCGStab | blockGMRES
    preconditioner  blockGAMG;       // blockDiagonal | blockGAMG
    tolerance       1e-8;            // absolute, on ||r||_2 / normFactor
    relTol          0.05;            // relative to initial residual of this outer iteration
    maxIter         200;
    gmresRestart    10;              // only blockGMRES
    blockGAMG
    {
        agglomerator      faceAreaPair;   // native GAMGAgglomeration
        nCellsInCoarsestLevel 500;
        mergeLevels       1;
        smoother          blockGaussSeidel;   // blockGaussSeidel | blockILU0
        nPreSweeps        1;
        nPostSweeps       2;
        nFinestSweeps     2;
        coarsestSolver    blockBiCGStab;
        coarsestTolerance 1e-3;
        maxOperatorComplexity 1.5;      // hard cap, see 6.3
        cacheAgglomeration true;
    }
}

Residual norm definition: ‖r‖ = sqrt( Σ_rows r_i² ) over all 4 components, with normFactor = Σ_i |A x_i| + |b_i| + SMALL computed once per outer iteration (analogue of native lduMatrix::solver::normFactor).

6.2 Block-BiCGStab / block-GMRES

Standard algorithms (Saad, Iterative Methods, Alg. 7.7 and 6.9) applied to the 4·nCells unknown vector, with the block matvec blockLduMatrix4::Amul() (LDU traversal, interface update as native lduMatrix::Amul with initMatrixInterfaces/updateMatrixInterfaces). Breakdown guards: if |ρ| < VSMALL·‖r₀‖² or |ω| < VSMALL → restart from current x (max 3 restarts, then return with converged=false).

6.3 Block-GAMG (additive correction, V-cycle)

Agglomeration: obtained from native GAMGAgglomeration::New(mesh, dict) on the p–p sub-block (faceAreaPair uses only geometry, so any sub-block works; use the p–p block's scalar coefficients if algebraicPair is selected). Reuse restrictAddressing, faceRestrictAddressing, nCells(level), and the native GAMGInterface objects for processor boundaries.

Coarse operators: Galerkin by summation (additive correction): coarse diag = Σ fine diag of agglomerated cells + Σ fine upper/lower of internal faces that become intra-coarse-cell; coarse upper/lower = Σ fine upper/lower of faces that survive. Same 4×4 layout. Sum in double, store as solveScalar.

Restriction: sum of residual blocks; prolongation: piecewise constant (add coarse correction to all fine cells of the coarse cell).

Smoother, block Gauss–Seidel: per cell, x_P ← D_P⁻¹ (b_P − Σ_N A_PN x_N), where D_P⁻¹ is the 4×4 dense inverse precomputed per level (LU with partial pivoting, block4Ops::invert, pivot guard: if |pivot| < 1e-30 add 1e-30 and count nSingularDiag for the log). Processor boundaries are treated Jacobi-style (values from the previous sweep), as native GaussSeidelSmoother.

Smoother, block ILU0 (optional): LDU-native incomplete block factorization with 4×4 block pivots; enable only if Gate B2 shows ≥ 15 % wall-clock gain on T2, otherwise keep as available option.

Operator complexity C_op = Σ_levels (nnz_blocks(level)) / nnz_blocks(fine). If C_op > maxOperatorComplexity, increase mergeLevels by 1 and re-agglomerate (max 3 attempts), then fail loudly. Log C_op, levels, cells per level once per run.

V-cycle as preconditioner (one V-cycle per Krylov iteration), not as a standalone solver.

6.4 Unit tests (Phase A gate, applications/test/)

Test-block4Ops: invert/LU of 1000 random well-conditioned 4×4 blocks; ‖A A⁻¹ − I‖_∞ < 1e-5 (float). Singular block → guard triggers, no NaN.

Test-blockMatrix: build a 4×4-block Poisson-like SPD system on the cavity mesh (diag = I·(Σ coeffs), off-diag = −I·coeff); compare Amul with the native lduMatrix::Amul applied component-wise: max abs diff < 1e-6·‖x‖. Run on 1 and 4 ranks; identical results to 1e-6.

Test-blockGAMG: solve the above system to 1e-8; iteration count ≤ 20; identical solution 1 vs 4 ranks to 1e-5.

Test-doubleReduce: sum of 10⁸ values of 1e-4 in float vs double-accumulated; the double result must equal 1e4 to 1e-9 relative; the float result must be visibly wrong (documents why the rule exists).

7. Adaptive convergence control

All quantities below are logged every iteration (Section 12).

7.1 Global CFL strategy (coupled.ptc)

Let R_n = ‖r_n‖ / ‖r_1‖ be the normalized initial residual of the coupled system at outer iteration n (before the solve), computed as in 6.1 (continuity row scaled).

cflStrategy mRDM (default):
CFL_{n+1} = min( CFL_max, CFL_n · max(1, min(β_max, (R_{n−1}/R_n)^γ)) ) if R_n ≤ R_{n−1}, else CFL_{n+1} = CFL_n (monotone: never decreases from residual stagnation). Defaults: γ = 1, β_max = 1.5.

cflStrategy EXP: CFL_{n+1} = min(CFL_max, CFL_n · β_exp), default β_exp = 1.1.

cflStrategy SER (available for comparison only): CFL_{n+1} = clamp(CFL_n · R_{n−1}/R_n, CFL_min, CFL_max).

Common: CFL_0 = 5, CFL_min = 1, CFL_max = 500. CFL is only decreased by the line search (7.2) or the sentinel (9.3), never by the strategy (except SER).

After any decrease, a hold of nHold = 10 iterations applies during which the strategy may not increase CFL.

7.2 Physicality line search (coupled.lineSearch)

Given increment dx = (dU, dp) from the linear solve:

ω_U = min_P ( 1, f_U · U_ref / max(|dU_P|, VSMALL) )

ω_p = min_P ( 1, f_p · p_ref / max(|dp_P|, VSMALL) )

ω = min(ω_U, ω_p), global minimum via returnReduce(..., minOp).

U_ref = max(|U_inlet|, gMax(|U|)) at iteration 1 and frozen; p_ref = ½ U_ref².

Defaults: f_U = 0.3, f_p = 0.5, omegaMin = 0.1, kappa = 0.5, maxCflCuts = 3 per outer iteration, beta = 1.0 (no additional CFL boost on full step; the strategy handles growth).

If ω < omegaMin: CFL ← CFL·κ, repeat the outer step (assemble again with the same fields but the new Δt). After maxCflCuts cuts: accept the step with ω = omegaMin and mark all cells with |dU_P| > f_U U_ref/omegaMin or |dp_P| > f_p p_ref/omegaMin into the dynamic remediation set (8.2).

7.3 Solution-limited local CFL (coupled.localLimit)

Estimated local change from the residual: δU_P = |r_P^{mom}| · Δt_P / V_P (r is the momentum residual of the current assembly). If δU_P > f_loc · U_ref (default f_loc = 0.5), cap Δt_P ← Δt_P · f_loc U_ref / δU_P for this cell only and record the count nLocalLimited. This keeps the global CFL aggressive while individual cells are throttled.

7.4 Start-up

potentialInit yes (Allrun calls native potentialFoam -writePhi -writep before coupledFoam; the solver reads the resulting U, p, phi).

startupUpwindIters 50: for the first 50 outer iterations β_P = 0 (pure upwind) everywhere; then β_P = 1 except remediated cells. The switch is also triggered earlier if R_n < 1e-2.

On restart, the start-up phase is skipped (state carries startupDone = true).

8. Remediation cell sets

Both sets are labelLists stored as cellSet objects (remediationStatic, remediationDynamic) and are written as a volScalarField remediationFlag (0 none, 1 static, 2 dynamic, 3 both) at every write time.

8.1 Static set (mesh quality, built once)

Cell is marked if any: non-orthogonality (max over its faces) > nonOrthThreshold (70°), skewness > skewThreshold (4), volume ratio to any neighbour > volRatioThreshold (20), aspect ratio > aspectThreshold (1000). Uses native primitiveMesh::checkFaceOrthogonality/checkFaceSkewness/… per-face results.
Treatment in static cells: β_P = β_static (default 0.0 = upwind), gradient limiter forced to cellLimited 1, non-orthogonal correction scaled by nonOrthLimiter_static (0.2), local CFL multiplied by cflFactor_static (0.5).

8.2 Dynamic set (solution spikes, every iteration)

Mark cell P if any of:

|U_P| > c_U · U_ref (default c_U = 2.0)

p_P < −c_p · p_ref or p_P > c_p · p_ref (default c_p = 5.0)

|U_P − Ū_N(P)| > c_spike · U_ref, with Ū_N(P) the arithmetic mean of face-neighbour values (default c_spike = 0.5)

any non-finite value in U_P, p_P, k_P, ω_P
Marked cells plus one layer of face neighbours (nLayers 1) receive: β_P = 0, local CFL factor cflFactor_dynamic (0.1), and — only if clipToNeighbourMean yes (default yes) — the increment dU_P is clipped so that |U_P + dU_P − Ū_N(P)| ≤ c_spike U_ref (increment clipping, never overwriting the field).
Hysteresis: a cell leaves the dynamic set after nQuietIters (20) consecutive iterations without re-marking. Warning: if size(dynamic ∪ static) > 0.01·nCells a warning is printed every 100 iterations; the run continues.

9. FPE safety and robustness

9.1 Build modes

Debug/Opt-dev runs (all tests in Section 13 except benchmarks): FOAM_SIGFPE=true FOAM_SETNAN=true. Foam::sigFpe enables FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW. Any trap = test failure.

Benchmark runs: FOAM_SIGFPE=false, plus _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON) and _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON) set at solver start when coupled.ftz yes (default yes in Opt). Log the mode.

9.2 Guard list (each must exist and be grep-able by the comment tag // GUARD:)

Every division by a computed quantity uses stabilise(x, SMALL) or an explicit max(x, VSMALL); grep / in src/ is reviewed at Gate A.

4×4 inversion pivot guard (6.3).

sqrt, log, pow arguments clamped to ≥ 0, ≥ VSMALL respectively.

k, ω floors after every correct() (5.8); ν_t capped at nutMax = 1e5·ν (log count if hit).

All reductions in double (6.0).

normFactor ≥ SMALL; ‖r₁‖ ≥ VSMALL before division in R_n.

Δt_P denominator (5.4) ≥ VSMALL.

Line-search denominators (7.2).

Row scaling s_p (5.3f) ≤ 1/VSMALL.

Conversion scalar → solveScalar at store time: values with |v| > 1e30 are clamped to ±1e30 and counted (nClampedCoeffs); a non-zero count is a warning in the log and a failure in tests T0–T3.

9.3 Sentinel and rollback

After each field update: one pass computing gMin/gMax of p, |U|, k, ω and nNonFinite (std::isfinite). Fail condition: nNonFinite > 0 or gMax(|U|) > 10 U_ref or |p| > 50 p_ref.
On fail: restore U, p, phi from the previous-iteration copies (which exist anyway for the PTC term; k, ω, ν_t are also copied — extra memory ≈ 3 scalar fields), set CFL ← max(CFL_min, 0.25·CFL), mark offending cells into the dynamic set, increment nRollbacks. After maxRollbacks (5) consecutive rollbacks: write the last valid fields to time <iter>_lastValid, write remediationFlag and a volScalarField sentinelFlag marking offending cells, print the offending cell centres (max 20), and exit with FatalError. Never write a time directory containing non-finite values.

10. Restart

On write: standard fields U, p, phi, k, omega, nut (binary, collated file handler recommended in Allrun via -fileHandler collated), plus coupledState (an IOdictionary in the time directory) containing: CFL, iter, startupDone, R_1 (initial residual reference), U_ref, p_ref, nHoldRemaining, dynamicSet (labelList), dynamicSetAge (labelList), refinementHistory (empty placeholder for Phase 2 AMR).

On start: if coupledState exists in startTime, all of the above are restored; else fresh start.

phi consistency: on read, phi is recomputed with 5.3(e) from the read U, p using the stored ā_P proxy (D_f reconstructed from the first assembly); the difference ‖phi_read − phi_recomputed‖_∞/‖phi‖_∞ is logged and must be < 1e-6 in test T-restart.

Atomic write: fields are written to <time>.tmp and renamed (native fileHandler writeObject then mv); implemented in coupledState::write().

11. Dictionary interface (system/fvSolution, subdict coupled) — complete keyword table

Keyword
Type
Default
Section

maxIter
label
2000
5.7

potentialInit
bool
yes
7.4

startupUpwindIters
label
50
7.4

nonOrthLimiter
scalar
0.5
5.3

ptc.cflStrategy
word
mRDM
7.1

ptc.CFL0 / CFLmin / CFLmax
scalar
5 / 1 / 500
7.1

ptc.gamma / betaMax / betaExp / nHold
scalar/label
1 / 1.5 / 1.1 / 10
7.1

lineSearch.fU / fp / omegaMin / kappa / maxCflCuts
scalar/label
0.3 / 0.5 / 0.1 / 0.5 / 3
7.2

localLimit.enabled / fLoc
bool/scalar
yes / 0.5
7.3

remediation.static.enabled
bool
yes
8.1

remediation.static.nonOrthThreshold / skewThreshold / volRatioThreshold / aspectThreshold
scalar
70 / 4 / 20 / 1000
8.1

remediation.static.beta / nonOrthLimiter / cflFactor
scalar
0 / 0.2 / 0.5
8.1

remediation.dynamic.enabled
bool
yes
8.2

remediation.dynamic.cU / cp / cSpike / nLayers / nQuietIters / cflFactor / clipToNeighbourMean
—
2 / 5 / 0.5 / 1 / 20 / 0.1 / yes
8.2

sentinel.maxRollbacks
label
5
9.3

bounds.kMin / omegaMin / nutMaxFactor
scalar
1e-12 / 1e-6 / 1e5
5.8, 9.2

convergence.forceCoeffsWindow / forceCoeffsTol / residualTol
label/scalar
100 / 0.002 / 1e-6
12.3

ftz
bool
yes
9.1

writeState
bool
yes
10

At start-up the solver prints the complete effective dictionary (defaults filled in) under the header coupledFoam: effective settings.

12. Logging, diagnostics, convergence

12.1 Per-iteration log line (single line, machine-parseable, prefix CF|)

CF| iter=<n> CFL=<v> omega=<v> cuts=<n> R=<Rn> rU=<v> rp=<v> linIters=<n> linRes=<v> tAsm=<s> tSolve=<s> tTurb=<s> tIter=<s> nStat=<n> nDyn=<n> nLocLim=<n> nRollback=<n> nClamped=<n>
rU, rp are the per-block normalized initial residuals; R the combined one used by 7.1.

12.2 Once-per-run block

Effective settings (11), mesh counts, static set size, GAMG levels/cells/C_op, precision check, FTZ state, peak RSS at end (from /proc/self/status VmHWM, gMax over ranks).

12.3 Convergence criterion (solver-internal stop)

Stop when both: (i) R_n < residualTol, or (ii) with a forceCoeffs function object present: over the last forceCoeffsWindow iterations, max(Cd) − min(Cd) ≤ forceCoeffsTol · |mean(Cd)| and the same for Cl. Criterion (ii) is the one used by the benchmark harness for both solvers (Section 14).

13. Test plan

All tests are run by pytest tests/ from the repository root; each test writes results/tests/<name>.json with the measured quantities and pass/fail. Native reference results are produced with simpleFoam on the identical mesh and schemes, converged to residualTol 1e-8 (SIMPLEC, consistent yes, tutorial relaxation factors). Unless stated, tests run on 4 ranks with scotch decomposition and on 1 rank; both must pass; single-precision cross-rank tolerance: 1e-4 relative on all integral quantities.

T0 — Lid-driven cavity (laminar), correctness of coupling

Mesh: tutorial incompressible/icoFoam/cavity/cavity blockMesh modified to 128×128×1. Re = 100 and 1000 (nu adjusted, U_lid = 1).

Solver setup: coupledFoam, laminar, schemes linearUpwind grad(U) for div, Gauss linear corrected for laplacian.

Pass: converge to R < 1e-8 within 300 iterations; centreline u(y) and v(x) profiles (sampled with sample at 129 points) vs. simpleFoam reference: L2 relative difference < 1e-3; no FPE trap in the Debug build; nClamped == 0.

T1 — pitzDaily (turbulent, SST) — PTC and turbulence coupling

Case: incompressible/simpleFoam/pitzDaily, turbulence switched to kOmegaSST, mesh as tutorial (~12 k cells).

Pass: converged (R < 1e-6) in ≤ 400 iterations; forceCoeffs-style integral check replaced by pressure drop inlet→outlet: within 1 % of simpleFoam reference; CFL ramp reached ≥ 100 without cuts after iteration 100; nRollbacks == 0.

T2 — backwardFacingStep2D (SST), reattachment

Case: incompressible/simpleFoam/backwardFacingStep2D (if absent in v2606, use pitzDailyExptInlet and skip the reattachment metric).

Pass: reattachment length (first sign change of wall shear on the lower wall downstream of the step, via wallShearStress FO) within 2 % of simpleFoam; ≥ 2× fewer outer iterations than simpleFoam to R < 1e-5.

T3 — airFoil2D (SST and GEKO) — external-aero-like, two models

Case: incompressible/simpleFoam/airFoil2D, run twice: kOmegaSST and GEKO (default coefficients).

Pass: Cl, Cd within 0.5 % of simpleFoam for each model; GEKO run completes with the same solver settings (no per-model tuning); nDyn at convergence == 0.

T-restart — exact resume (run on T1 and T3-SST)

Run 150 iterations, write; restart from 150 to convergence. Compare with the uninterrupted run: iteration count to convergence identical ±2; final Cd/Cl/Δp identical to 1e-5 relative; phi-consistency metric (Section 10) < 1e-6.

T-fpe — torture test (run on T1)

Init U = 0 everywhere, potentialInit no, CFL0 200, startupUpwindIters 0.

Pass in Debug build with traps: no trap; run recovers (nRollbacks ≥ 1 allowed) and converges to R < 1e-5. Pass in Opt build: identical converged Δp to 1e-4.

T4 — motorBike (snappyHexMesh, real bad cells, MRF none)

Case: incompressible/simpleFoam/motorBike, two meshes: (a) tutorial refinement (~350 k cells), (b) refined: snappyHexMeshDict refinement levels +1 on the motorBike surface and refinement box, target 3–5 M cells; report actual count. 16 ranks.

Pass: converges by criterion 12.3(ii) with forceCoeffsWindow 100, tol 0.002; Cd within 1 % of simpleFoam on the same mesh; static set ≤ 1 % of cells; nRollbacks == 0; peak RSS recorded.

T5 — Ahmed body 25° slant (generated geometry)

Geometry: cases/T5_ahmed/geometry/make_ahmed.py generates ahmed25.stl (numpy-stl): length 1044 mm, width 389 mm, height 288 mm, front edges rounded R = 100 mm (both vertical and horizontal front edges), rear slant length 222 mm at 25°, ground clearance 50 mm, no stilts. Domain: inlet 3 L upstream, outlet 7 L downstream, half-width 2 W (symmetry plane at y = 0, half model), height 5 H. Ground: moving-wall fixedValue U = U_inf. U_inf = 40 m/s, ν = 1.5e-5. Turbulence: SST, nutkWallFunction, y⁺ target 30–100 (layers: 5, expansion 1.2).

Mesh: snappyHexMesh, target 5–8 M cells (half model); report actual.

Pass: converged by 12.3(ii); Cd within 1 % of simpleFoam on the same mesh; both within ±10 % of the experimental value 0.285 (Ahmed et al. 1984) — the experimental comparison is informational, the solver-to-solver one is the pass criterion; nRollbacks == 0; peak RSS recorded.

T-scaling — strong scaling (on T4b)

Ranks 1, 2, 4, 8, 16, 300 iterations each, both solvers. Record time/iteration. Pass: parallel efficiency of coupledFoam at 16 ranks ≥ 80 % of simpleFoam's efficiency at 16 ranks (i.e. the coupled solver must not scale worse than native).

14. Benchmark protocol (bench/run_bench.py)

For each case in {T1, T2, T3-SST, T3-GEKO, T4a, T4b, T5}:

Identical mesh, decomposition (scotch), schemes (linearUpwind grad(U), Gauss linear corrected), turbulence model and BCs for all solvers.

Solvers compared: (A) simpleFoam with tutorial settings (SIMPLE, p: GAMG, U/k/ω: smoothSolver); (B) simpleFoam SIMPLEC (consistent yes, relaxation p 1.0 / U 0.9 / k,ω 0.9) — the tuned native baseline; (C) coupledFoam defaults; (D) coupledFoam with preconditioner blockDiagonal (isolates the AMG benefit).

Each configuration is run 3 times; report median and min/max. Every rank is wrapped by bench/rank_wrapper.sh (/usr/bin/time -v) → per-rank max RSS, summed and max.

Stop criterion for all solvers: 12.3(ii) with window 100 and tol 0.002 evaluated by the harness on the forceCoeffs (or Δp for T1/T2) history — the solvers' own criteria are disabled (maxIter large) so the harness decides identically.

Metrics per run: iters_to_conv, wall_to_conv_s, time_per_iter_s, t_assembly, t_linsolve, t_turb (from CF| lines; for simpleFoam from -profiling or wall/iter only), peakRSS_GB_sum, peakRSS_GB_max_rank, Cd_final, Cl_final, speedup = wall_B / wall_C.

Output: results/bench/<case>_<solver>_<run>.json and a consolidated results/bench/summary.csv.

15. Report specification (bench/make_report.py → report/REPORT.md, English)

Sections, in this order, all generated from results/:

Executive summary — one table: case, cells, ranks, native-best wall-clock, coupled wall-clock, speedup, memory ratio, Cd agreement.

Method — 1 page: block coupling, PTC, block-GAMG, remediation, precision (link to this spec).

Correctness — T0 profile plots (coupled vs native), T3 Cl/Cd table, T5 Cd vs experiment.

Convergence behaviour — per case: residual vs iteration (both solvers), Cd vs iteration with the convergence window marked, CFL history and ω history of the coupled solver.

Performance — bar charts of wall-clock to convergence (A/B/C/D), time/iteration breakdown (assembly / linear solve / turbulence), strong-scaling plot (T-scaling), memory bar chart with the Section 1 budget line.

Robustness — T-fpe outcome, remediation set sizes over iterations for T4/T5, rollback counts.

Linear solver — block-GAMG levels, C_op, Krylov iterations per outer iteration; blockDiagonal vs blockGAMG comparison.

Limitations and Phase 2 — explicit list from DECISIONS.md.

Reproduction — exact commands and commit hash.
Figures: PNG, 150 dpi, report/figures/, matplotlib, consistent colours (native = grey, coupled = blue).

16. Implementation phases and gates

Phase
Content
Gate (all mandatory)

0
Toolchain: two OpenFOAM v2606 SPDP builds (Opt, Debug), Test-precision, repo skeleton, CI script Allrun-tests.sh
sizeof(solveScalar)==4; simpleFoam pitzDaily tutorial runs in both builds; pytest tests/test_env.py green

A
blockMatrix, block4Ops, doubleReduce, block interfaces, blockBiCGStab, blockDiagonal precon; unit tests 6.4
Test-block4Ops, Test-blockMatrix, Test-doubleReduce green on 1 and 4 ranks; guard review (9.2 item 1) recorded in results/gates/phase_A.json

B1
Assembly 5.3 (no MRF), BCs 5.5 (fixedValue/zeroGradient/slip/outlet), flux update, outer loop 5.7 with fixed CFL, sentinel 9.3, turbulence 5.8
T0 green (both Re), T1 green with fixed CFL 20; Debug-build no traps

B2
blockGAMG (GS smoother), C_op cap, optional ILU0 with decision
Test-blockGAMG green; T1 wall-clock with blockGAMG < 0.5× blockDiagonal; ILU0 decision in DECISIONS.md

C
PTC 5.4, CFL strategies 7.1, line search 7.2, local limit 7.3, start-up 7.4
T1, T2 green incl. CFL criteria; T-fpe green

D
Remediation 8, restart 10, coupledState, MRF 5.6, inletOutlet/totalPressure BCs
T3 (SST+GEKO), T-restart green; T4a green

E
T4b, T5 geometry + meshing + runs, T-scaling, benchmark harness, report
T4b, T5, T-scaling green; report/REPORT.md generated with all figures; every number in the executive summary traceable to a JSON

Estimated effort ordering (largest first): B2, B1, C, D, E, A. If time-boxed, E is never cut — a smaller T5 mesh is acceptable (minimum 2 M cells) but the report must exist.

17. Phase 2 (not part of the one-shot; design constraints already honoured)

SAMG (Uroić/Jasak): selective coarsening on the p–p block, ILUC0 smoother; plug-in as another blockSolver preconditioner.

AMR: dynamicRefineFvMesh, hex-only, protection layers, refinement every 200–500 iterations, fvMeshDistribute rebalancing, memory guard; coupledState.refinementHistory already reserved.

petsc4Foam FieldSplit reference (single-precision PETSc).

Implicit cyclic/AMI coupling in the block matrix.

45 M-cell production run: user acceptance step with the memory budget of v1 (100–115 GB) verified via peak RSS.

18. References

Darwish M., Sraj I., Moukalled F., A coupled incompressible flow solver on structured grids, Numer. Heat Transfer B 52 (2007).

Mangani L., Buchmayr M., Darwish M., Development of a Novel Fully Coupled Solver in OpenFOAM: Steady-State Incompressible Turbulent Flows, Numer. Heat Transfer B 66 (2014); … in Rotational Reference Frames (2014); A fully coupled OpenFOAM solver for transient incompressible turbulent flows in ALE formulation (2017) — implicit boundary treatment.

Uroić T., Jasak H., Block-selective algebraic multigrid for implicitly coupled pressure-velocity system, Comput. Fluids (2018); Parallelisation of selective algebraic multigrid for block–pressure–velocity system in OpenFOAM, Comput. Phys. Commun. (2020).

Oliani S., Carnevale M. et al., ICSFoam: An OpenFOAM library for implicit coupled simulations of high-speed flows, Comput. Phys. Commun. (2023); GPU follow-up arXiv:2403.07882 (2024).

Ceze M., Fidkowski K., A Robust Adaptive Solution Strategy for High-Order Implicit CFD Solvers, AIAA 2011-3696; Pseudo-transient Continuation, Solution Update Methods, and CFL Strategies for DG Discretizations of the RANS-SA Equations, AIAA 2013 — SER/EXP/RDM/mRDM, physicality line search.

Kelley C.T., Keyes D.E., Convergence analysis of pseudo-transient continuation, SIAM J. Numer. Anal. (1998).

Solution-limited time stepping to enhance reliability in CFD applications, J. Comput. Phys. 228 (2009).

Clain S., Diot S., Loubère R., A high-order finite volume method for systems of conservation laws — Multi-dimensional Optimal Order Detection (MOOD), J. Comput. Phys. (2011).

Saad Y., Iterative Methods for Sparse Linear Systems, 2nd ed., SIAM (2003).

Ahmed S.R., Ramm G., Faltin G., Some salient features of the time-averaged ground vehicle wake, SAE 840300 (1984).

OpenFOAM v2606 release notes: GEKO (Menter & Matyushenko 2025), AMI cache restart.

