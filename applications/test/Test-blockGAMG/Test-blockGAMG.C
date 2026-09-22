/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-blockGAMG

Description
    Unit test of the block-GAMG-preconditioned block Krylov solver
    (spec 6.4, amendment B9): solve the 4x4-block Poisson-like system of
    testPoissonSystem.H with b = A x_exact to tolerance 1e-8 (or tighter).

    Pass: converged within 20 iterations. The 1-vs-N-rank comparison (1e-5,
    done by pytest) uses global integrals of the solution:
    sum_P x_k V_P and sum_P x_k^2 V_P per component, and the max error to
    x_exact, or the full solution written with -dumpSolution.

    The same solve with the blockDiagonal preconditioner is run for
    comparison (informational, skipped with -skipDiagonal).

    Solver settings: solvers.coupled of system/fvSolution if present,
    otherwise the spec 6.1 defaults; tolerance 1e-8 (or -tolerance) and
    relTol 0 always.

    B9 additions:
      -cycle V|F|W|K     sets blockGAMG.cycleType (default: the dictionary
                         value, else K) and switches autoTune off so that
                         the requested cycle is the one measured.
                         cycleType K forces solver blockFGMRES (B1).
      -mergeLevels n     sets blockGAMG.mergeLevels
      -procAgglom on|off sets blockGAMG.processorAgglomerator to
                         masterCoarsest / none (6.3.4, Test-procAgglom)
    The JSON records cycleType, mergeLevels (requested and used), the
    measured coarsening ratios, nIterations, wall seconds, converged.

Usage
    Test-blockGAMG [-parallel] [-json <file>]
        [-solver blockBiCGStab|blockGMRES|blockFGMRES] [-tolerance <tol>]
        [-cycle V|F|W|K] [-mergeLevels <n>] [-procAgglom on|off]
        [-dumpSolution <file>] [-skipDiagonal]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "testPoissonSystem.H"
#include "blockSolver.H"
#include "blockGAMGPrecon.H"
#include "doubleReduce.H"
#include "jsonWriter.H"
#include "clockTime.H"

using namespace Foam;

// Test parameters (fixed by spec 6.4)
static constexpr doubleScalar testTolerance = 1e-8;
static constexpr label maxAllowedIters = 20;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct gamgInfo
{
    bool valid = false;
    word cycleType;
    label mergeLevels = -1;
    label nLevels = 0;
    doubleScalar Cop = 0;
    List<doubleScalar> ratios;
};


struct solveResult
{
    blockScalarList x;
    blockSolverPerformance perf;
    FixedList<doubleScalar, 4> intX = FixedList<doubleScalar, 4>(Zero);
    FixedList<doubleScalar, 4> intX2 = FixedList<doubleScalar, 4>(Zero);
    doubleScalar maxErr = 0;
    doubleScalar seconds = 0;

    //- ||b - A x||/normFactor of the returned x, evaluated independently of
    //  the solver (informational: a solver may report a recursively updated
    //  residual instead of the true one)
    doubleScalar trueResidual = 0;
};


static solveResult runSolve
(
    const fvMesh& mesh,
    const blockLduMatrix4& B,
    const blockScalarList& xExact,
    const dictionary& solverDict,
    dictionary& settings,
    dictionary& gamgStats,
    gamgInfo& info
)
{
    blockScalarList b(B.nRows());
    B.Amul(b, xExact);

    blockScalarList x(B.nRows(), Zero);

    autoPtr<blockSolver> solver = blockSolver::New(B, solverDict);
    solver->writeSettings(settings);

    clockTime timer;
    solveResult res;
    res.perf = solver->solve(x, b);
    res.seconds = timer.elapsedTime();

    const blockGAMGPrecon* gp =
        dynamic_cast<const blockGAMGPrecon*>(solver->preconditioner());
    if (gp)
    {
        const blockGAMG& g = gp->gamg();
        gamgStats = g.statsDict();
        info.valid = true;
        info.cycleType = blockGAMG::cycleName(g.cycleType());
        info.mergeLevels = g.mergeLevels();
        info.nLevels = g.nLevels();
        info.Cop = g.operatorComplexity();
        const List<reduceScalar>& r = g.ratios();
        info.ratios.resize(r.size());
        forAll(r, i)
        {
            info.ratios[i] = r[i];
        }
    }

    res.x = x;

    // True residual of the returned x, accumulated in double, with the
    // solver's normFactor for the zero initial guess (sum |A 0| + |b| + SMALL)
    {
        List<reduceScalar> xd(x.size());
        forAll(x, i)
        {
            xd[i] = toDouble(x[i]);
        }
        blockScalarList r(B.nRows());
        B.residualDouble(r, xd, b);
        const doubleScalar nf =
            doubleReduce::sumMag(b, B.comm()) + doubleScalarSMALL;
        res.trueResidual = doubleReduce::norm2(r, B.comm())/nf;
    }

    const scalarField& V = mesh.V();
    res.intX = Zero;
    res.intX2 = Zero;
    res.maxErr = 0;
    forAll(V, celli)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            const doubleScalar xv = toDouble(x[celli*blockDim + k]);
            res.intX[k] += xv*V[celli];
            res.intX2[k] += xv*xv*V[celli];
            res.maxErr = std::max
            (
                res.maxErr,
                std::abs(xv - toDouble(xExact[celli*blockDim + k]))
            );
        }
    }
    doubleReduce::parSum(res.intX.data(), 4, UPstream::worldComm);
    doubleReduce::parSum(res.intX2.data(), 4, UPstream::worldComm);
    Foam::reduce(res.maxErr, maxOp<doubleScalar>());

    return res;
}


int main(int argc, char *argv[])
{
    argList::addOption("json", "file", "Write results as JSON");
    argList::addOption
    (
        "solver",
        "name",
        "blockBiCGStab (default) | blockGMRES | blockFGMRES;"
        " cycleType K forces blockFGMRES"
    );
    argList::addOption
    (
        "tolerance",
        "value",
        "Solver tolerance (default 1e-8, spec 6.4; tighter values allowed)"
    );
    argList::addOption
    (
        "dumpSolution",
        "file",
        "Write the blockGAMG solution gathered on the master as CSV"
        " (x y z x0 x1 x2 x3), for the 1-vs-N-rank comparison"
    );
    argList::addOption
    (
        "cycle",
        "V|F|W|K",
        "blockGAMG cycleType (default: dictionary value, else K);"
        " also sets autoTune off"
    );
    argList::addOption("mergeLevels", "n", "blockGAMG mergeLevels");
    argList::addOption
    (
        "procAgglom",
        "on|off",
        "on: processorAgglomerator masterCoarsest (6.3.4);"
        " off: no processor agglomeration (default: dictionary)"
    );
    argList::addBoolOption
    (
        "skipDiagonal",
        "Skip the informational blockDiagonal comparison solve"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    tmp<volScalarField> tpsi = testPoisson::makePsi(mesh);
    tmp<fvScalarMatrix> tM = testPoisson::makeScalarMatrix(tpsi.ref());

    blockLduMatrix4 B(mesh);
    testPoisson::fillBlockMatrix(tM(), B);

    const blockScalarList xExact = testPoisson::makeX(mesh);

    // Solver controls
    dictionary solverDict(mesh.solversDict().subOrEmptyDict("coupled"));

    const word requestedSolver =
        args.getOrDefault<word>("solver", "blockBiCGStab");
    solverDict.set("solver", requestedSolver);
    solverDict.set("preconditioner", word("blockGAMG"));
    const doubleScalar tol =
        args.getOrDefault<doubleScalar>("tolerance", testTolerance);
    if (tol > testTolerance)
    {
        FatalErrorInFunction
            << "tolerance " << tol << " looser than the spec value "
            << testTolerance << exit(FatalError);
    }
    solverDict.set("tolerance", tol);
    solverDict.set("relTol", doubleScalar(0));

    // blockGAMG sub-dictionary: cycle type, mergeLevels, proc agglomeration
    dictionary gamgDict(solverDict.subOrEmptyDict("blockGAMG"));

    word cycleType = gamgDict.getOrDefault<word>("cycleType", "K");
    if (args.found("cycle"))
    {
        cycleType = args.get<word>("cycle");
        // The measured cycle must be the requested one (6.3.5)
        gamgDict.set("autoTune", false);
    }
    if
    (
        cycleType != "V" && cycleType != "F"
     && cycleType != "W" && cycleType != "K"
    )
    {
        FatalErrorInFunction
            << "cycleType " << cycleType << " is not one of V F W K"
            << exit(FatalError);
    }
    gamgDict.set("cycleType", cycleType);

    const label mergeLevelsRequested =
        args.getOrDefault<label>
        (
            "mergeLevels",
            gamgDict.getOrDefault<label>("mergeLevels", 1)
        );
    if (args.found("mergeLevels"))
    {
        gamgDict.set("mergeLevels", mergeLevelsRequested);
    }

    word procAgglom("dictionary");
    if (args.found("procAgglom"))
    {
        procAgglom = args.get<word>("procAgglom");
        if (procAgglom == "on")
        {
            gamgDict.set("processorAgglomerator", word("masterCoarsest"));
        }
        else if (procAgglom == "off")
        {
            gamgDict.set("processorAgglomerator", word("none"));
        }
        else
        {
            FatalErrorInFunction
                << "-procAgglom " << procAgglom << ": expected on or off"
                << exit(FatalError);
        }
    }
    if (requestedSolver != "blockFGMRES")
    {
        // The unit test compares the V/F/W cycles with the requested
        // (fixed-preconditioner) solver as before; the iterative coarsest
        // solve is accepted explicitly here (D-069 F10)
        gamgDict.set("allowVariableCoarsest", true);
    }
    solverDict.set("blockGAMG", gamgDict);

    // B1: cycleType K requires blockFGMRES (variable preconditioner)
    bool solverForced = false;
    if (cycleType == "K" && requestedSolver != "blockFGMRES")
    {
        solverDict.set("solver", word("blockFGMRES"));
        solverForced = true;
        Info<< "cycleType K: solver " << requestedSolver
            << " replaced by blockFGMRES (amendment B1)" << nl;
    }

    Info<< "blockGAMG test: cycleType " << cycleType
        << ", mergeLevels " << mergeLevelsRequested
        << ", procAgglom " << procAgglom
        << ", solver " << solverDict.get<word>("solver") << nl << endl;

    dictionary settingsGAMG, statsGAMG;
    gamgInfo infoG;
    const solveResult rG =
        runSolve(mesh, B, xExact, solverDict, settingsGAMG, statsGAMG, infoG);

    Info<< "blockGAMG:     ";
    rG.perf.print(Info);

    // Informational: blockDiagonal with the requested (not forced) solver
    const bool runDiagonal = !args.found("skipDiagonal");
    solveResult rD;
    if (runDiagonal)
    {
        dictionary diagDict(solverDict);
        diagDict.set("solver", requestedSolver);
        diagDict.set("preconditioner", word("blockDiagonal"));
        diagDict.set("maxIter", label(10000));
        dictionary settingsDiag, statsDiag;
        gamgInfo infoD;
        rD = runSolve
        (
            mesh, B, xExact, diagDict, settingsDiag, statsDiag, infoD
        );

        Info<< "blockDiagonal: ";
        rD.perf.print(Info);
    }

    const bool pass =
        rG.perf.converged && rG.perf.nIterations <= maxAllowedIters;

    Info<< "GAMG stats " << statsGAMG << nl
        << "cycleType  " << infoG.cycleType << nl
        << "mergeLevels used " << infoG.mergeLevels << nl
        << "ratios     " << infoG.ratios << nl
        << "int x      " << rG.intX << nl
        << "int x^2    " << rG.intX2 << nl
        << "max |x - x_exact| " << rG.maxErr << nl
        << "true residual of the returned x " << rG.trueResidual << nl
        << "wall GAMG " << rG.seconds << " s, diagonal " << rD.seconds
        << " s" << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("dumpSolution"))
    {
        // Gather cell centres and solution values on the master
        List<pointField> allC(UPstream::nProcs());
        List<List<blockScalar>> allX(UPstream::nProcs());
        allC[UPstream::myProcNo()] = mesh.C().primitiveField();
        allX[UPstream::myProcNo()] = rG.x;
        Pstream::gatherList(allC);
        Pstream::gatherList(allX);

        if (UPstream::master())
        {
            const fileName f(args.get<fileName>("dumpSolution"));
            mkDir(f.path());
            OFstream os(f);
            os.precision(12);
            forAll(allC, proci)
            {
                const pointField& C = allC[proci];
                const List<blockScalar>& X = allX[proci];
                forAll(C, i)
                {
                    os  << C[i].x() << ' ' << C[i].y() << ' ' << C[i].z();
                    for (label k = 0; k < blockDim; ++k)
                    {
                        os  << ' ' << toDouble(X[i*blockDim + k]);
                    }
                    os  << nl;
                }
            }
        }
    }

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-blockGAMG");
        j.add("nProcs", UPstream::nProcs());
        j.add("nCells", returnReduce(mesh.nCells(), sumOp<label>()));
        j.add("solver", solverDict.get<word>("solver"));
        j.add("solverRequested", requestedSolver);
        j.add("solverForcedByCycleK", solverForced);
        j.add("tolerance", tol);
        j.add("cycleType", infoG.valid ? infoG.cycleType : cycleType);
        j.add("cycleTypeRequested", cycleType);
        j.add("mergeLevelsRequested", mergeLevelsRequested);
        j.add("mergeLevels", infoG.mergeLevels);
        j.add("procAgglom", procAgglom);
        j.addList("ratios", infoG.ratios);
        j.add("converged", rG.perf.converged);
        j.add("nIterations", rG.perf.nIterations);
        j.add("maxAllowedIterations", maxAllowedIters);
        j.add("initialResidual", rG.perf.initialResidual);
        j.add("finalResidual", rG.perf.finalResidual);
        j.add("rho", rG.perf.rho);
        j.addList("intX", List<doubleScalar>(rG.intX));
        j.addList("intX2", List<doubleScalar>(rG.intX2));
        j.add("maxErrorToExact", rG.maxErr);
        j.add("wallSeconds", rG.seconds);
        j.add("gamgLevels", infoG.nLevels);
        j.add("gamgCop", infoG.Cop);
        j.addList
        (
            "gamgCellsPerLevel",
            statsGAMG.getOrDefault<labelList>("cellsPerLevel", labelList())
        );
        j.add("diagonalRun", runDiagonal);
        j.add("diagonalConverged", rD.perf.converged);
        j.add("diagonalIterations", rD.perf.nIterations);
        j.add("diagonalWallSeconds", rD.seconds);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
