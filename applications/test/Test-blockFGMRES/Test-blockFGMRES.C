/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-blockFGMRES

Description
    Unit test of the flexible GMRES (amendment B9): with a deliberately
    variable preconditioner FGMRES must converge; standard (right-
    preconditioned) GMRES on the same set-up is allowed to stall. This
    documents why the K-cycle (a variable preconditioner) requires
    blockFGMRES (B1).

    System: the 4x4-block Poisson-like system of testPoissonSystem.H on the
    case mesh (T0 cavity), b = A x_exact, x0 = 0.

    Variable preconditioner "testVariableGS" (defined and registered in the
    runtime selection table of blockPreconditioner by this executable):
    every application computes w = M_j^-1 r as n_j forward block
    Gauss-Seidel sweeps (Foam::blockGaussSeidel) from w = 0, with n_j drawn
    uniformly from nMin..nMax by a Foam::Random with a fixed seed (the same
    sequence on every rank, so the sweeps stay collective in parallel).
    Standard GMRES forms x += M^-1 (V y) once per restart cycle with a new
    n_j, i.e. with an operator that did not generate the Krylov basis.

    Runs (same restart, maxIter, tolerance, relTol 0, seed):
      fgmres       blockFGMRES + variable preconditioner   (must converge)
      gmres        blockGMRES  + variable preconditioner   (recorded only)
      gmresFixed   blockGMRES  + fixed nMax sweeps          (control: shows
                   that the set-up itself is solvable by GMRES)

    For every run the true residual ||b - A x||_2 / (sum |b| + SMALL) (the
    solver norm of spec 6.1 for x0 = 0) is recomputed after the solve and
    compared with the residual the solver reported.

    Pass: fgmres converged and its true residual <= trueResidualFactor
    (10) * tolerance.

Usage
    Test-blockFGMRES [-parallel] [-json <file>] [-tolerance 1e-8]
        [-restart 50] [-maxIter 5000] [-nMin 1] [-nMax 6] [-seed 1234]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "testPoissonSystem.H"
#include "blockSolver.H"
#include "blockPreconditioner.H"
#include "blockGaussSeidel.H"
#include "doubleReduce.H"
#include "jsonWriter.H"
#include "clockTime.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class testVariableGS Declaration
\*---------------------------------------------------------------------------*/

//- Variable preconditioner: a random number of block Gauss-Seidel sweeps
//  per application. Settings from the solver dict's "testVariableGS"
//  sub-dictionary: nMin, nMax, seed.
class testVariableGS
:
    public blockPreconditioner
{
    // Private Data

        const dictionary coeffs_;

        const label nMin_;
        const label nMax_;
        const label seed_;

        //- The smoother
        blockGaussSeidel gs_;

        //- Random sweep counts (identical sequence on every rank)
        mutable Random rnd_;

        //- Statistics
        mutable label nApplications_;
        mutable label nSweepsTotal_;
        mutable label nSweepsMinUsed_;
        mutable label nSweepsMaxUsed_;


public:

    //- Runtime type information
    TypeName("testVariableGS");


    // Constructors

        testVariableGS(const blockSolver& solver, const dictionary& dict)
        :
            blockPreconditioner(solver),
            coeffs_(dict.subOrEmptyDict("testVariableGS")),
            nMin_(max(label(1), coeffs_.getOrDefault<label>("nMin", 1))),
            nMax_(max(nMin_, coeffs_.getOrDefault<label>("nMax", 6))),
            seed_(coeffs_.getOrDefault<label>("seed", 1234)),
            gs_(solver.matrix(), coeffs_),
            rnd_(seed_),
            nApplications_(0),
            nSweepsTotal_(0),
            nSweepsMinUsed_(labelMax),
            nSweepsMaxUsed_(0)
        {}


    //- Destructor
    virtual ~testVariableGS() = default;


    // Member Functions

        label nApplications() const noexcept { return nApplications_; }
        label nSweepsTotal() const noexcept { return nSweepsTotal_; }
        label nSweepsMinUsed() const noexcept { return nSweepsMinUsed_; }
        label nSweepsMaxUsed() const noexcept { return nSweepsMaxUsed_; }

        virtual void update()
        {
            gs_.update();
        }

        virtual void precondition
        (
            blockScalarUList& w,
            const blockScalarUList& r
        ) const
        {
            label n = nMin_;
            if (nMax_ > nMin_)
            {
                const doubleScalar u = rnd_.sample01<scalar>();
                n = nMin_
                  + static_cast<label>
                    (
                        std::floor(u*doubleScalar(nMax_ - nMin_ + 1))
                    );
                n = min(max(n, nMin_), nMax_);
            }

            w = Zero;
            gs_.smooth(w, r, n);

            ++nApplications_;
            nSweepsTotal_ += n;
            nSweepsMinUsed_ = min(nSweepsMinUsed_, n);
            nSweepsMaxUsed_ = max(nSweepsMaxUsed_, n);
        }

        virtual label nSingularDiag() const
        {
            return gs_.nSingularDiag();
        }

        virtual void writeSettings(dictionary& dict) const
        {
            dict.add("type", type());
            dict.add("nMin", nMin_);
            dict.add("nMax", nMax_);
            dict.add("seed", seed_);
        }
};


defineTypeNameAndDebug(testVariableGS, 0);
addToRunTimeSelectionTable(blockPreconditioner, testVariableGS, dictionary);

} // End namespace Foam


using namespace Foam;

// Pass criterion: true residual within this factor of the tolerance
static constexpr doubleScalar trueResidualFactor = 10;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct runResult
{
    word name;
    word solver;
    bool variable = true;
    blockSolverPerformance perf;
    doubleScalar trueResidual = 0;
    doubleScalar maxErr = 0;
    doubleScalar seconds = 0;
    label nApplications = 0;
    doubleScalar meanSweeps = 0;
    label minSweeps = 0;
    label maxSweeps = 0;
};


static runResult runSolve
(
    const word& name,
    const blockLduMatrix4& B,
    const blockScalarList& b,
    const blockScalarList& xExact,
    const dictionary& solverDict,
    const bool variable
)
{
    runResult res;
    res.name = name;
    res.solver = solverDict.get<word>("solver");
    res.variable = variable;

    blockScalarList x(B.nRows(), Zero);

    autoPtr<blockSolver> solver = blockSolver::New(B, solverDict);

    clockTime timer;
    res.perf = solver->solve(x, b);
    res.seconds = timer.elapsedTime();

    const testVariableGS* pc =
        dynamic_cast<const testVariableGS*>(solver->preconditioner());
    if (pc)
    {
        res.nApplications = pc->nApplications();
        res.meanSweeps =
            pc->nApplications() > 0
          ? doubleScalar(pc->nSweepsTotal())/doubleScalar(pc->nApplications())
          : 0;
        res.minSweeps = pc->nApplications() > 0 ? pc->nSweepsMinUsed() : 0;
        res.maxSweeps = pc->nSweepsMaxUsed();
    }

    // True residual in the solver norm (x0 = 0: normFactor = sum |b|)
    blockScalarList Ax(B.nRows());
    B.Amul(Ax, x);
    blockScalarList r(B.nRows());
    forAll(r, i)
    {
        r[i] = b[i] - Ax[i];
    }
    FixedList<reduceScalar, 2> s;
    s[0] = doubleReduce::localSumSqr(r);
    s[1] = doubleReduce::localSumMag(b);
    doubleReduce::parSum(s.data(), 2, B.comm());
    res.trueResidual = std::sqrt(s[0])/(s[1] + doubleScalarSMALL);

    forAll(x, i)
    {
        res.maxErr = std::max
        (
            res.maxErr,
            std::abs(toDouble(x[i]) - toDouble(xExact[i]))
        );
    }
    Foam::reduce(res.maxErr, maxOp<doubleScalar>());

    Info<< name << ": ";
    res.perf.print(Info);
    Info<< "    true residual " << res.trueResidual
        << ", max |x - x_exact| " << res.maxErr
        << ", preconditioner applications " << res.nApplications
        << " (sweeps " << res.minSweeps << ".." << res.maxSweeps
        << ", mean " << res.meanSweeps << ")"
        << ", wall " << res.seconds << " s" << nl << endl;

    return res;
}


static void addRun(jsonWriter& j, const runResult& r)
{
    jsonWriter o;
    o.add("solver", r.solver);
    o.add("variablePreconditioner", r.variable);
    o.add("converged", r.perf.converged);
    o.add("breakdown", r.perf.breakdown);
    o.add("nIterations", r.perf.nIterations);
    o.add("nRestarts", r.perf.nRestarts);
    o.add("initialResidual", r.perf.initialResidual);
    o.add("finalResidual", r.perf.finalResidual);
    o.add("trueResidual", r.trueResidual);
    o.add("maxErrorToExact", r.maxErr);
    o.add("preconditionerApplications", r.nApplications);
    o.add("sweepsMin", r.minSweeps);
    o.add("sweepsMax", r.maxSweeps);
    o.add("sweepsMean", r.meanSweeps);
    o.add("wallSeconds", r.seconds);
    j.addRaw(r.name, o.str());
}


int main(int argc, char *argv[])
{
    argList::addOption("json", "file", "Write results as JSON");
    argList::addOption("tolerance", "value", "Solver tolerance (1e-8)");
    argList::addOption("restart", "m", "(F)GMRES restart (50)");
    argList::addOption("maxIter", "n", "Maximum iterations (5000)");
    argList::addOption("nMin", "n", "Minimum sweeps per application (1)");
    argList::addOption("nMax", "n", "Maximum sweeps per application (6)");
    argList::addOption("seed", "n", "Random seed (1234)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    tmp<volScalarField> tpsi = testPoisson::makePsi(mesh);
    tmp<fvScalarMatrix> tM = testPoisson::makeScalarMatrix(tpsi.ref());

    blockLduMatrix4 B(mesh);
    testPoisson::fillBlockMatrix(tM(), B);

    const blockScalarList xExact = testPoisson::makeX(mesh);
    blockScalarList b(B.nRows());
    B.Amul(b, xExact);

    const doubleScalar tol = args.getOrDefault<doubleScalar>("tolerance", 1e-8);
    const label restart = args.getOrDefault<label>("restart", 50);
    const label maxIter = args.getOrDefault<label>("maxIter", 5000);
    const label nMin = args.getOrDefault<label>("nMin", 1);
    const label nMax = args.getOrDefault<label>("nMax", 6);
    const label seed = args.getOrDefault<label>("seed", 1234);

    // Common controls (no Eisenstat-Walker, relTol 0)
    dictionary base;
    base.set("tolerance", tol);
    base.set("relTol", doubleScalar(0));
    base.set("adaptiveRelTol", false);
    base.set("minIter", label(1));
    base.set("maxIter", maxIter);
    base.set("restart", restart);
    base.set("preconditioner", word("testVariableGS"));
    {
        dictionary c;
        c.set("nMin", nMin);
        c.set("nMax", nMax);
        c.set("seed", seed);
        base.set("testVariableGS", c);
    }

    dictionary fDict(base);
    fDict.set("solver", word("blockFGMRES"));

    dictionary gDict(base);
    gDict.set("solver", word("blockGMRES"));

    dictionary gFixedDict(gDict);
    {
        dictionary c;
        c.set("nMin", nMax);
        c.set("nMax", nMax);
        c.set("seed", seed);
        gFixedDict.set("testVariableGS", c);
    }

    Info<< "Test-blockFGMRES: tolerance " << tol << ", restart " << restart
        << ", maxIter " << maxIter << ", sweeps " << nMin << ".." << nMax
        << ", seed " << seed << nl << endl;

    const runResult rF = runSolve("fgmres", B, b, xExact, fDict, true);
    const runResult rG = runSolve("gmres", B, b, xExact, gDict, true);
    const runResult rGF =
        runSolve("gmresFixed", B, b, xExact, gFixedDict, false);

    const bool fgmresOk =
        rF.perf.converged && rF.trueResidual <= trueResidualFactor*tol;
    const bool gmresTrulyConverged =
        rG.perf.converged && rG.trueResidual <= trueResidualFactor*tol;
    const bool pass = fgmresOk;

    word gmresOutcome;
    if (gmresTrulyConverged)
    {
        gmresOutcome = "converged";
    }
    else if (rG.perf.converged)
    {
        gmresOutcome = "falseConvergence";
    }
    else
    {
        gmresOutcome = "stalled";
    }

    Info<< "FGMRES (variable):     "
        << (fgmresOk ? "converged" : "NOT converged") << nl
        << "GMRES (variable):      " << gmresOutcome << nl
        << "GMRES (fixed control): "
        << (rGF.perf.converged ? "converged" : "NOT converged") << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-blockFGMRES");
        j.add("nProcs", UPstream::nProcs());
        j.add("nCells", returnReduce(mesh.nCells(), sumOp<label>()));
        j.add("tolerance", tol);
        j.add("restart", restart);
        j.add("maxIter", maxIter);
        j.add("sweepsMin", nMin);
        j.add("sweepsMax", nMax);
        j.add("seed", seed);
        j.add("trueResidualFactor", trueResidualFactor);
        addRun(j, rF);
        addRun(j, rG);
        addRun(j, rGF);
        j.add("fgmresConverged", fgmresOk);
        j.add("gmresOutcome", gmresOutcome);
        j.add("gmresFixedConverged", rGF.perf.converged);
        j.add
        (
            "note",
            "Standard GMRES forms x += M^-1 (V y) with a preconditioner that"
            " did not generate the Krylov basis; with a variable M this is"
            " inconsistent, so GMRES may stall or report false convergence."
            " FGMRES stores Z_j = M_j^-1 v_j and stays exact (amendment B1)."
        );
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
