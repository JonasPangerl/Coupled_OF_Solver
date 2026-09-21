/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-blockGAMG

Description
    Unit test of the block-GAMG-preconditioned block Krylov solver
    (spec 6.4): solve the 4x4-block Poisson-like system of
    testPoissonSystem.H with b = A x_exact to tolerance 1e-8.

    Pass: converged within 20 iterations. The 1-vs-N-rank comparison (1e-5,
    done by pytest) uses global integrals of the solution:
    sum_P x_k V_P and sum_P x_k^2 V_P per component, and the max error to
    x_exact.

    The same solve with the blockDiagonal preconditioner is run for
    comparison (informational).

    Solver settings: solvers.coupled of system/fvSolution if present,
    otherwise the spec 6.1 defaults with tolerance 1e-8 and relTol 0.

Usage
    Test-blockGAMG [-parallel] [-json <file>] [-solver blockBiCGStab|blockGMRES]

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

struct solveResult
{
    blockSolverPerformance perf;
    FixedList<doubleScalar, 4> intX;
    FixedList<doubleScalar, 4> intX2;
    doubleScalar maxErr;
    doubleScalar seconds;
};


static solveResult runSolve
(
    const fvMesh& mesh,
    const blockLduMatrix4& B,
    const blockScalarList& xExact,
    const dictionary& solverDict,
    dictionary& settings,
    dictionary& gamgStats
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
        gamgStats = gp->gamg().statsDict();
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
    argList::addOption("solver", "name", "blockBiCGStab (default) | blockGMRES");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    tmp<volScalarField> tpsi = testPoisson::makePsi(mesh);
    tmp<fvScalarMatrix> tM = testPoisson::makeScalarMatrix(tpsi.ref());

    blockLduMatrix4 B(mesh);
    testPoisson::fillBlockMatrix(tM(), B);

    const blockScalarList xExact = testPoisson::makeX(mesh);

    // Solver controls
    dictionary solverDict;
    {
        const dictionary& fvSol = mesh.solverDict("coupled");
        solverDict = fvSol;
    }
    solverDict.set("solver", args.getOrDefault<word>("solver", "blockBiCGStab"));
    solverDict.set("preconditioner", word("blockGAMG"));
    solverDict.set("tolerance", testTolerance);
    solverDict.set("relTol", doubleScalar(0));

    dictionary settingsGAMG, statsGAMG;
    const solveResult rG =
        runSolve(mesh, B, xExact, solverDict, settingsGAMG, statsGAMG);

    Info<< "blockGAMG:     ";
    rG.perf.print(Info);

    dictionary diagDict(solverDict);
    diagDict.set("preconditioner", word("blockDiagonal"));
    diagDict.set("maxIter", label(10000));
    dictionary settingsDiag, statsDiag;
    const solveResult rD =
        runSolve(mesh, B, xExact, diagDict, settingsDiag, statsDiag);

    Info<< "blockDiagonal: ";
    rD.perf.print(Info);

    const bool pass =
        rG.perf.converged && rG.perf.nIterations <= maxAllowedIters;

    Info<< "GAMG stats " << statsGAMG << nl
        << "int x      " << rG.intX << nl
        << "int x^2    " << rG.intX2 << nl
        << "max |x - x_exact| " << rG.maxErr << nl
        << "wall GAMG " << rG.seconds << " s, diagonal " << rD.seconds
        << " s" << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-blockGAMG");
        j.add("nProcs", UPstream::nProcs());
        j.add("nCells", returnReduce(mesh.nCells(), sumOp<label>()));
        j.add("solver", solverDict.get<word>("solver"));
        j.add("converged", rG.perf.converged);
        j.add("nIterations", rG.perf.nIterations);
        j.add("maxAllowedIterations", maxAllowedIters);
        j.add("initialResidual", rG.perf.initialResidual);
        j.add("finalResidual", rG.perf.finalResidual);
        j.addList("intX", List<doubleScalar>(rG.intX));
        j.addList("intX2", List<doubleScalar>(rG.intX2));
        j.add("maxErrorToExact", rG.maxErr);
        j.add("wallSeconds", rG.seconds);
        j.add("gamgLevels", statsGAMG.getOrDefault<label>("nLevels", 0));
        j.add("gamgCop", statsGAMG.getOrDefault<doubleScalar>("Cop", 0));
        j.addList
        (
            "gamgCellsPerLevel",
            statsGAMG.getOrDefault<labelList>("cellsPerLevel", labelList())
        );
        j.add("diagonalConverged", rD.perf.converged);
        j.add("diagonalIterations", rD.perf.nIterations);
        j.add("diagonalWallSeconds", rD.seconds);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
