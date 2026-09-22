/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockSolver.H"
#include "blockPreconditioner.H"
#include "doubleReduce.H"
#include "blockKernels.H"
#include "coupledDefaults.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockSolver, 0);
    defineRunTimeSelectionTable(blockSolver, dictionary);
}


// * * * * * * * * * * * * blockSolverPerformance  * * * * * * * * * * * * * //

void Foam::blockSolverPerformance::print(Ostream& os) const
{
    os  << solverName
        << ": initial residual = " << initialResidual
        << ", final residual = " << finalResidual
        << ", no. iterations " << nIterations;
    if (nRestarts)
    {
        os  << ", restarts " << nRestarts;
    }
    if (breakdown)
    {
        os  << ", BREAKDOWN";
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockSolver::blockSolver
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    matrix_(matrix),
    controlDict_(dict),
    tolerance_
    (
        dict.getOrDefault<doubleScalar>("tolerance", coupledDefaults::tolerance)
    ),
    relTol_(dict.getOrDefault<doubleScalar>("relTol", coupledDefaults::relTol)),
    maxIter_(dict.getOrDefault<label>("maxIter", coupledDefaults::maxIter)),
    minIter_(dict.getOrDefault<label>("minIter", coupledDefaults::minIter)),
    maxRestarts_
    (
        dict.getOrDefault<label>("maxRestarts", coupledDefaults::maxRestarts)
    ),
    preconPtr_(),
    preconVersion_(-1),
    diag_(nullptr)
{
    const word preconName =
        dict.getOrDefault<word>("preconditioner", "none");

    if (preconName != "none")
    {
        preconPtr_ = blockPreconditioner::New(*this, dict);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blockSolver> Foam::blockSolver::New
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
{
    const word solverType(dict.get<word>("solver"));

    auto* ctorPtr = dictionaryConstructorTable(solverType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "blockSolver",
            solverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<blockSolver>(ctorPtr(matrix, dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockSolver::~blockSolver()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::blockSolver::converged
(
    const reduceScalar residual,
    const reduceScalar initialResidual,
    const label nIter
) const
{
    if (nIter < minIter_)
    {
        return false;
    }
    return
    (
        residual < tolerance_
     || (relTol_ > 0 && residual < relTol_*initialResidual)
    );
}


void Foam::blockSolver::updatePreconditioner() const
{
    if (preconPtr_ && preconVersion_ != matrix_.version())
    {
        if (diagActive(1))
        {
            const doubleScalar t0 = diagnostics::clock();
            preconPtr_->update();
            diag_->addPrecSetup(diagnostics::clock() - t0);
        }
        else
        {
            preconPtr_->update();
        }
        preconVersion_ = matrix_.version();
    }
}


void Foam::blockSolver::precondition
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    if (preconPtr_)
    {
        if (diagActive(1))
        {
            if (diag_->active(2))
            {
                diag_->beginApplication();
            }
            const doubleScalar t0 = diagnostics::clock();
            preconPtr_->precondition(w, r);
            diag_->addPrecApply(diagnostics::clock() - t0);
            if (diag_->active(2))
            {
                diag_->endApplication();
            }
        }
        else
        {
            preconPtr_->precondition(w, r);
        }
    }
    else
    {
        blockKernels::copy(r.size(), r.cdata(), w.data());
    }
}


void Foam::blockSolver::diagBeginSolve() const
{
    if (diagActive(2))
    {
        diag_->beginSolve();
    }
}


void Foam::blockSolver::diagEndSolve(const blockSolverPerformance& perf) const
{
    if (diagActive(2))
    {
        diag_->endSolve
        (
            perf.solverName.c_str(),
            perf.initialResidual,
            perf.finalResidual,
            perf.nIterations,
            perf.nRestarts,
            perf.converged,
            perf.rho
        );
    }
}


Foam::reduceScalar Foam::blockSolver::normFactor
(
    const blockScalarUList& x,
    const blockScalarUList& b
) const
{
    blockScalarList Ax(matrix_.nRows());
    matrix_.Amul(Ax, x);

    reduceScalar s = 0;
    const label n = matrix_.nRows();
    const blockScalar* __restrict__ AxPtr = Ax.cdata();
    const blockScalar* __restrict__ bPtr = b.cdata();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += std::abs(toDouble(AxPtr[i])) + std::abs(toDouble(bPtr[i]));
    }

    // GUARD: normFactor >= SMALL (spec 9.2)
    return doubleReduce::parSum(s, matrix_.comm()) + doubleScalarSMALL;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::reduceScalar Foam::blockSolver::measureRho
(
    const blockScalarUList& r
) const
{
    updatePreconditioner();
    const label n = matrix_.nRows();
    blockScalarList z(n), Az(n);
    precondition(z, r);
    matrix_.Amul(Az, z);
    for (label i = 0; i < n; ++i)
    {
        Az[i] = r[i] - Az[i];
    }
    // GUARD: ||r|| > 0
    return
        doubleReduce::norm2(Az, matrix_.comm())
       /max(doubleReduce::norm2(r, matrix_.comm()), doubleScalarVSMALL);
}


void Foam::blockSolver::writeSettings(dictionary& dict) const
{
    dict.add("solver", type());
    dict.add("tolerance", tolerance_);
    dict.add("relTol", relTol_);
    dict.add("maxIter", maxIter_);
    dict.add("minIter", minIter_);
    dict.add("maxRestarts", maxRestarts_);
    if (preconPtr_)
    {
        dictionary pd;
        preconPtr_->writeSettings(pd);
        dict.add("preconditioner", pd);
    }
    else
    {
        dict.add("preconditioner", word("none"));
    }
}


// ************************************************************************* //
