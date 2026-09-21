# F1 half-car readiness check (HANDOFF 4.6, D-031)

Reviewed 2026-09-21 (read-only) against the user's case repository,
`f1_halfcar/runs/s5_amr` (the most recent stage; solver simpleFoam,
kOmegaSST, ~10 ranks). Feature-by-feature support status of coupledFoam:

| F1 case feature | coupledFoam status |
|---|---|
| inlet: `fixedValue` U, `zeroGradient` p | supported (generic BC linearisation, D-002; exercised in T1-T3) |
| outlet: `inletOutlet` U, `fixedValue` p | supported; the flux-sign switch enters via valueInternalCoeffs and is lagged one outer iteration, same as simpleFoam |
| `slip` tunnelWall (U) | supported through the native diagonal linearisation (the tangential-projection off-diagonals are lagged, exactly as in simpleFoam). On the review list to double-check |
| `noSlip` car/driver, `fixedValue` moving road | supported |
| `rotatingWallVelocity` wheels | supported (value-type BC, generic path) |
| symmetry plane, constraint patches via `setConstraintTypes` | symmetry/symmetryPlane/wedge are dedicated zero-flux patches; processor patches are implicit block interfaces |
| MRF: 2 wheel zones, rear with `nonRotatingPatches` | MRFCoupling (D-006). `nonRotatingPatches` handling is on the review list - verify before the first F1 run |
| kOmegaSST | supported: segregated turbulence after each coupled step (as T3/T4) |
| `div(phi,U) bounded Gauss linearUpwind limited` | supported: the momentum operator is built from the native fvVectorMatrix; the deferred correction blends upwind -> case scheme with beta |
| no cyclicAMI, no in-solver mesh changes (AMR is a meshing stage) | OK. coupledFoam does not support runtime topology change |
| open domain (outlet fixes p) | good: the closed-domain pressure-reference path (D-021) is not needed |
| function objects `forcesPerPatch`, `solverInfo`, post-processing | forces work (T3-T5 use forceCoeffs). `solverInfo` sees no native per-field solver data from the block solve; expect empty entries there (the CF\| log line and the summary JSON carry the solver data) |

Conclusion: nothing in the case is outside coupledFoam's supported set.
Two review items to close before the first production F1 run: the `slip`
linearisation and MRF `nonRotatingPatches` (both assigned to the assembly
review).

## Drop-in changes for the F1 case

1. `system/controlDict`: `application coupledFoam;`
2. `system/fvSolution`: replace the `solvers` p/U entries and the
   `SIMPLE`/`relaxationFactors` blocks by:

```
solvers
{
    coupled
    {
        solver          blockFGMRES;
        preconditioner  blockGAMG;
        tolerance       1e-8;
        adaptiveRelTol  yes;            // Eisenstat-Walker (B2)
        minIter         1;
        maxIter         200;
        restart         10;
        blockGAMG
        {
            agglomerator          faceAreaPair;
            mergeLevels           2;
            nCellsInCoarsestLevel 200;
            cycleType             K;    // use V until the K-cycle fix lands
            smoother              blockILU0;  // pending gate B2 on T2
            nPreSweeps            1;
            nPostSweeps           2;
            nFinestSweeps         2;
            coarsestSolver        blockBiCGStab;
            coarsestTolerance     1e-3;
            coarsestMaxIter       50;
            maxOperatorComplexity 1.5;
            cacheAgglomeration    true;
            autoTune              yes;
            // processorAgglomerator masterCoarsest; (default, D-030)
        }
    }

    // k/omega stay segregated: keep the existing smoothSolver entries.
}

coupled
{
    maxIter         2000;
    potentialInit   yes;
    ptc
    {
        // defaults (coupledDefaults.H); CFLmax etc. tune later
    }
    convergence
    {
        residualTol     1e-6;   // start value; tighten per force history
        forceCoeffsWindow 100;
        forceCoeffsTol  0.002;
    }
}
```

3. Initialise with `potentialFoam -writePhi -writep` (the `potentialInit
   yes` above expects that; our cases/*/Allrun show the call order).
4. Keep decomposePar/scotch as-is; run
   `mpirun --bind-to core --map-by core -np 10 coupledFoam -parallel`
   (D-031: at most 10 cores on this machine).

Before recommending this for production: finish HANDOFF 4.1 (parallel
MPI fix + K cycle), the T2 smoother gate (B2), and T3/T4 validation, then
re-check this file.
