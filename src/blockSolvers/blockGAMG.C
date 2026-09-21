/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGAMG.H"
#include "doubleReduce.H"
#include "blockPreconditioner.H"
#include "coupledDefaults.H"
#include "lduMesh.H"
#include "objectRegistry.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGAMG::blockGAMG
(
    const blockLduMatrix4& fine,
    const dictionary& dict
)
:
    fine_(fine),
    dict_(dict),
    nPreSweeps_
    (
        dict.getOrDefault<label>("nPreSweeps", coupledDefaults::nPreSweeps)
    ),
    nPostSweeps_
    (
        dict.getOrDefault<label>("nPostSweeps", coupledDefaults::nPostSweeps)
    ),
    nFinestSweeps_
    (
        dict.getOrDefault<label>
        (
            "nFinestSweeps",
            coupledDefaults::nFinestSweeps
        )
    ),
    maxCop_
    (
        dict.getOrDefault<doubleScalar>
        (
            "maxOperatorComplexity",
            coupledDefaults::maxOperatorComplexity
        )
    ),
    maxCopAttempts_
    (
        dict.getOrDefault<label>
        (
            "maxCopAttempts",
            coupledDefaults::maxCopAttempts
        )
    ),
    aggPtr_(nullptr),
    mergeLevelsUsed_(-1),
    Cop_(0),
    coarse_(),
    smoothers_(),
    coarsestSolver_(),
    x_(),
    b_(),
    r_(),
    nCoarsestIters_(0)
{
    // Effective defaults written back so that the effective settings print
    // shows every value actually used
    dict_.add
    (
        "agglomerator",
        dict.getOrDefault<word>("agglomerator", "faceAreaPair"),
        false
    );
    dict_.add
    (
        "nCellsInCoarsestLevel",
        dict.getOrDefault<label>
        (
            "nCellsInCoarsestLevel",
            coupledDefaults::nCellsInCoarsestLevel
        ),
        false
    );
    dict_.add
    (
        "mergeLevels",
        dict.getOrDefault<label>("mergeLevels", coupledDefaults::mergeLevels),
        false
    );
    dict_.add
    (
        "smoother",
        dict.getOrDefault<word>("smoother", "blockGaussSeidel"),
        false
    );
    dict_.add
    (
        "coarsestSolver",
        dict.getOrDefault<word>("coarsestSolver", "blockBiCGStab"),
        false
    );
    dict_.add
    (
        "coarsestTolerance",
        dict.getOrDefault<doubleScalar>
        (
            "coarsestTolerance",
            coupledDefaults::coarsestTolerance
        ),
        false
    );
    dict_.add
    (
        "coarsestMaxIter",
        dict.getOrDefault<label>
        (
            "coarsestMaxIter",
            coupledDefaults::coarsestMaxIter
        ),
        false
    );
    dict_.add
    (
        "cacheAgglomeration",
        dict.getOrDefault<bool>
        (
            "cacheAgglomeration",
            coupledDefaults::cacheAgglomeration
        ),
        false
    );

    agglomerate();

    const GAMGAgglomeration& agg = *aggPtr_;
    const label nCoarse = agg.size();

    // Coarse matrices on the native coarse lduMesh levels
    coarse_.resize(nCoarse);
    for (label l = 0; l < nCoarse; ++l)
    {
        const lduMesh& cMesh = agg.meshLevel(l + 1);
        coarse_.set
        (
            l,
            new blockLduMatrix4
            (
                cMesh.lduAddr(),
                agg.interfaceLevel(l + 1),
                cMesh.comm()
            )
        );
    }

    // Smoothers on all but the coarsest level
    const label nLev = nLevels();
    smoothers_.resize(nLev - 1);
    for (label l = 0; l < nLev - 1; ++l)
    {
        smoothers_.set(l, blockSmoother::New(matrixLevel(l), dict_).ptr());
    }

    // Coarsest-level solver
    {
        dictionary cd;
        cd.add("solver", dict_.get<word>("coarsestSolver"));
        cd.add("preconditioner", word("blockDiagonal"));
        cd.add("tolerance", doubleScalar(0));
        cd.add("relTol", dict_.get<doubleScalar>("coarsestTolerance"));
        cd.add("maxIter", dict_.get<label>("coarsestMaxIter"));
        cd.add
        (
            "pivotGuard",
            dict_.getOrDefault<doubleScalar>
            (
                "pivotGuard",
                coupledDefaults::pivotGuard
            )
        );
        coarsestSolver_ = blockSolver::New(matrixLevel(nLev - 1), cd);
    }

    // Work vectors (level 0 uses the caller's vectors)
    x_.resize(nLev);
    b_.resize(nLev);
    r_.resize(nLev);
    for (label l = 0; l < nLev; ++l)
    {
        const label n = matrixLevel(l).nRows();
        if (l > 0)
        {
            x_[l].resize(n, Zero);
            b_[l].resize(n, Zero);
        }
        r_[l].resize(n, Zero);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::reduceScalar Foam::blockGAMG::operatorComplexity
(
    const GAMGAgglomeration& agg
) const
{
    reduceScalar nnz[2] = {reduceScalar(fine_.nnzBlocks()), 0};
    nnz[1] = nnz[0];

    for (label l = 0; l < agg.size(); ++l)
    {
        reduceScalar n = reduceScalar(agg.nCells(l)) + 2*reduceScalar(agg.nFaces(l));

        const lduInterfacePtrsList coarseIfaces = agg.interfaceLevel(l + 1);
        const labelList& nPatchFaces = agg.nPatchFaces(l);
        forAll(coarseIfaces, i)
        {
            if
            (
                coarseIfaces.set(i)
             && blockLduInterface::isBlockCoupled(coarseIfaces[i])
            )
            {
                n += reduceScalar(nPatchFaces[i]);
            }
        }
        nnz[1] += n;
    }

    doubleReduce::parSum(nnz, 2, fine_.comm());

    // GUARD: a mesh has at least one cell globally
    return nnz[1]/std::max(nnz[0], reduceScalar(1));
}


void Foam::blockGAMG::agglomerate()
{
    const lduMesh& mesh = fine_.mesh();
    label m = dict_.get<label>("mergeLevels");

    for (label attempt = 0; attempt < maxCopAttempts_; ++attempt, ++m)
    {
        dictionary aggDict(dict_);
        aggDict.set("mergeLevels", m);
        aggDict.set("name", word("blockGAMGAgglomeration_m" + Foam::name(m)));

        const GAMGAgglomeration& agg = GAMGAgglomeration::New(mesh, aggDict);

        const reduceScalar cop = operatorComplexity(agg);

        Info<< "blockGAMG: agglomeration attempt " << attempt + 1
            << " mergeLevels " << m << " levels " << agg.size() + 1
            << " C_op " << cop << endl;

        if (cop <= maxCop_)
        {
            aggPtr_ = &agg;
            mergeLevelsUsed_ = m;
            Cop_ = cop;
            return;
        }
    }

    FatalErrorInFunction
        << "Operator complexity exceeds maxOperatorComplexity " << maxCop_
        << " after " << maxCopAttempts_ << " agglomeration attempts"
        << " (last mergeLevels " << m - 1 << ")." << nl
        << "Increase mergeLevels or nCellsInCoarsestLevel (spec 6.3)."
        << exit(FatalError);
}


void Foam::blockGAMG::restrictMatrix(const label fineLevel)
{
    const GAMGAgglomeration& agg = *aggPtr_;
    const blockLduMatrix4& F = matrixLevel(fineLevel);
    blockLduMatrix4& C = coarse_[fineLevel];

    const label nCC = agg.nCells(fineLevel);
    const label nCF = agg.nFaces(fineLevel);

    // Galerkin sums in double (spec 6.3)
    reduceScalarList cDiag(blockSize*nCC, Zero);
    reduceScalarList cUpper(blockSize*nCF, Zero);
    reduceScalarList cLower(blockSize*nCF, Zero);

    const labelField& restrictAddr = agg.restrictAddressing(fineLevel);
    const labelList& faceRestrictAddr = agg.faceRestrictAddressing(fineLevel);
    const boolList& flip = agg.faceFlipMap(fineLevel);
    const bool haveFlip = (flip.size() == faceRestrictAddr.size());

    const blockScalarList& fD = F.diag();
    const blockScalarList& fU = F.upper();
    const blockScalarList& fL = F.lower();

    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockSize; ++k)
        {
            cDiag[c*blockSize + k] += toDouble(fD[celli*blockSize + k]);
        }
    }

    forAll(faceRestrictAddr, facei)
    {
        const label cFace = faceRestrictAddr[facei];

        if (cFace >= 0)
        {
            const bool flipped = haveFlip && flip[facei];
            const blockScalar* fu = fU.cdata() + facei*blockSize;
            const blockScalar* fl = fL.cdata() + facei*blockSize;
            reduceScalar* cu = cUpper.data() + cFace*blockSize;
            reduceScalar* cl = cLower.data() + cFace*blockSize;

            for (label k = 0; k < blockSize; ++k)
            {
                if (!flipped)
                {
                    cu[k] += toDouble(fu[k]);
                    cl[k] += toDouble(fl[k]);
                }
                else
                {
                    cu[k] += toDouble(fl[k]);
                    cl[k] += toDouble(fu[k]);
                }
            }
        }
        else
        {
            // Face internal to a coarse cell: both couplings go to the
            // coarse diagonal (the cell couples to itself)
            reduceScalar* cd = cDiag.data() + (-1 - cFace)*blockSize;
            const blockScalar* fu = fU.cdata() + facei*blockSize;
            const blockScalar* fl = fL.cdata() + facei*blockSize;
            for (label k = 0; k < blockSize; ++k)
            {
                cd[k] += toDouble(fu[k]) + toDouble(fl[k]);
            }
        }
    }

    // Store narrowed
    blockScalarList& cD = C.diag();
    blockScalarList& cU = C.upper();
    blockScalarList& cL = C.lower();
    forAll(cD, i)
    {
        cD[i] = narrow(cDiag[i]);
    }
    forAll(cU, i)
    {
        cU[i] = narrow(cUpper[i]);
        cL[i] = narrow(cLower[i]);
    }

    // Processor interface coefficients
    const labelListList& patchFaceRestrict =
        agg.patchFaceRestrictAddressing(fineLevel);

    forAll(C.interfaces(), ci)
    {
        const label inti = C.interfaces()[ci].index();

        // Matching fine block interface (same native interface index)
        label fi = -1;
        forAll(F.interfaces(), j)
        {
            if (F.interfaces()[j].index() == inti)
            {
                fi = j;
                break;
            }
        }
        if (fi < 0)
        {
            FatalErrorInFunction
                << "No fine-level block interface for coarse interface "
                << inti << " on level " << fineLevel + 1
                << abort(FatalError);
        }

        const blockScalarList& fC = F.interfaceCoeffs(fi);
        blockScalarList& cC = C.interfaceCoeffs(ci);
        const labelList& pfr = patchFaceRestrict[inti];

        reduceScalarList sumC(cC.size(), Zero);
        forAll(pfr, pf)
        {
            const label cpf = pfr[pf];
            for (label k = 0; k < blockSize; ++k)
            {
                sumC[cpf*blockSize + k] += toDouble(fC[pf*blockSize + k]);
            }
        }
        forAll(cC, i)
        {
            cC[i] = narrow(sumC[i]);
        }
    }

    C.markUpdated();
}


void Foam::blockGAMG::restrictResidual(const label l) const
{
    const labelField& restrictAddr = aggPtr_->restrictAddressing(l);
    blockScalarList& bc = b_[l + 1];
    const blockScalarList& rf = r_[l];

    bc = Zero;
    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockDim; ++k)
        {
            bc[c*blockDim + k] += rf[celli*blockDim + k];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::blockGAMG::globalCellsPerLevel() const
{
    labelList n(nLevels());
    n[0] = fine_.nCells();
    for (label l = 1; l < nLevels(); ++l)
    {
        n[l] = coarse_[l - 1].nCells();
    }
    Foam::reduce(n.data(), int(n.size()), sumOp<label>(), UPstream::msgType(), fine_.comm());
    return n;
}


Foam::label Foam::blockGAMG::nSingularDiag() const
{
    label n = 0;
    forAll(smoothers_, l)
    {
        n += smoothers_[l].nSingularDiag();
    }
    if (coarsestSolver_ && coarsestSolver_->preconditioner())
    {
        n += coarsestSolver_->preconditioner()->nSingularDiag();
    }
    return n;
}


void Foam::blockGAMG::update()
{
    if (!dict_.get<bool>("cacheAgglomeration"))
    {
        // Re-agglomerate: drop the cached object and select again. The level
        // sizes must not change for a static mesh; guard anyway.
        const label nLevOld = nLevels();
        fine_.mesh().thisDb().checkOut
        (
            const_cast<GAMGAgglomeration*>(aggPtr_)
        );
        aggPtr_ = nullptr;
        agglomerate();
        if (aggPtr_->size() + 1 != nLevOld)
        {
            FatalErrorInFunction
                << "Level count changed on re-agglomeration of a static mesh"
                << abort(FatalError);
        }
    }

    for (label l = 0; l < coarse_.size(); ++l)
    {
        restrictMatrix(l);
    }

    forAll(smoothers_, l)
    {
        smoothers_[l].update();
    }
}


void Foam::blockGAMG::Vcycle
(
    blockScalarUList& x,
    const blockScalarUList& b
) const
{
    const label nLev = nLevels();

    auto X = [&](const label l) -> blockScalarUList&
    {
        return (l == 0 ? x : static_cast<blockScalarUList&>(x_[l]));
    };
    auto B = [&](const label l) -> const blockScalarUList&
    {
        return (l == 0 ? b : static_cast<const blockScalarUList&>(b_[l]));
    };

    // Downward leg
    for (label l = 0; l < nLev - 1; ++l)
    {
        X(l) = Zero;
        if (nPreSweeps_ > 0)
        {
            smoothers_[l].smooth(X(l), B(l), nPreSweeps_);
        }
        matrixLevel(l).residual(r_[l], X(l), B(l));
        restrictResidual(l);
    }

    // Coarsest level
    {
        blockScalarUList& xc = X(nLev - 1);
        xc = Zero;
        const blockSolverPerformance perf =
            coarsestSolver_->solve(xc, B(nLev - 1));
        nCoarsestIters_ = perf.nIterations;
    }

    // Upward leg
    for (label l = nLev - 2; l >= 0; --l)
    {
        // Prolong: X(l) += P x_{l+1}
        const labelField& restrictAddr = aggPtr_->restrictAddressing(l);
        const blockScalarList& xc = x_[l + 1];
        blockScalarUList& xf = X(l);
        forAll(restrictAddr, celli)
        {
            const label c = restrictAddr[celli];
            for (label k = 0; k < blockDim; ++k)
            {
                xf[celli*blockDim + k] += xc[c*blockDim + k];
            }
        }

        smoothers_[l].smooth
        (
            xf,
            B(l),
            (l == 0 ? nFinestSweeps_ : nPostSweeps_)
        );
    }
}


void Foam::blockGAMG::writeStats(Ostream& os) const
{
    const labelList n = globalCellsPerLevel();
    os  << "blockGAMG: levels " << nLevels()
        << ", mergeLevels " << mergeLevelsUsed_
        << ", C_op " << Cop_ << nl
        << "blockGAMG: cells per level " << flatOutput(n) << endl;
}


Foam::dictionary Foam::blockGAMG::statsDict() const
{
    dictionary d;
    d.add("nLevels", nLevels());
    d.add("mergeLevels", mergeLevelsUsed_);
    d.add("Cop", Cop_);
    d.add("cellsPerLevel", globalCellsPerLevel());
    return d;
}


void Foam::blockGAMG::writeSettings(dictionary& dict) const
{
    dict.merge(dict_);
    dict.set("nPreSweeps", nPreSweeps_);
    dict.set("nPostSweeps", nPostSweeps_);
    dict.set("nFinestSweeps", nFinestSweeps_);
    dict.set("maxOperatorComplexity", maxCop_);
    dict.set("maxCopAttempts", maxCopAttempts_);
}


// ************************************************************************* //
