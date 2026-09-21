/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGAMG.H"
#include "blockPreconditioner.H"
#include "block4Ops.H"
#include "doubleReduce.H"
#include "coupledDefaults.H"
#include "lduMesh.H"
#include "objectRegistry.H"
#include "PstreamReduceOps.H"
#include <cmath>

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::blockGAMG::cycleKind Foam::blockGAMG::cycleFromWord(const word& w)
{
    if (w == "V") return cycleKind::V;
    if (w == "F") return cycleKind::F;
    if (w == "W") return cycleKind::W;
    if (w == "K") return cycleKind::K;
    FatalErrorInFunction
        << "Unknown cycleType " << w << ", valid: V F W K"
        << exit(FatalError);
    return cycleKind::V;
}


Foam::word Foam::blockGAMG::cycleName(const cycleKind c)
{
    switch (c)
    {
        case cycleKind::V: return "V";
        case cycleKind::F: return "F";
        case cycleKind::W: return "W";
        case cycleKind::K: return "K";
    }
    return "?";
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGAMG::blockGAMG
(
    const blockLduMatrix4& fine,
    const dictionary& dict
)
:
    fine_(fine),
    dict_(dict),
    cycle_
    (
        cycleFromWord(dict.getOrDefault<word>("cycleType", "K"))
    ),
    kThreshold_
    (
        dict.getOrDefault<doubleScalar>
        (
            "kCycleThreshold",
            coupledDefaults::kCycleThreshold
        )
    ),
    kMaxSteps_
    (
        dict.getOrDefault<label>("kCycleMaxSteps", coupledDefaults::kCycleMaxSteps)
    ),
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
    minRatio_
    (
        dict.getOrDefault<doubleScalar>
        (
            "minCoarseningRatio",
            coupledDefaults::minCoarseningRatio
        )
    ),
    maxMergeLevels_
    (
        dict.getOrDefault<label>
        (
            "maxMergeLevels",
            coupledDefaults::maxMergeLevels
        )
    ),
    denseLUMaxCells_
    (
        dict.getOrDefault<label>
        (
            "denseLUMaxCells",
            coupledDefaults::denseLUMaxCells
        )
    ),
    aggPtr_(nullptr),
    mergeLevelsUsed_(-1),
    Cop_(0),
    ratios_(),
    coarse_(),
    smoothers_(),
    coarsestSolver_(),
    useDenseLU_(false),
    denseLU_(),
    densePivot_(),
    nCoarsestIters_(0)
{
    // Effective defaults written back so that the effective-settings print
    // shows every value actually used
    auto setDefault = [this, &dict](const word& key, const auto& value)
    {
        if (!dict.found(key))
        {
            dict_.add(key, value);
        }
    };
    setDefault("agglomerator", word("faceAreaPair"));
    setDefault("nCellsInCoarsestLevel", coupledDefaults::nCellsInCoarsestLevel);
    setDefault("mergeLevels", coupledDefaults::mergeLevels);
    setDefault("smoother", word("blockGaussSeidel"));
    setDefault("coarsestSolver", word("blockBiCGStab"));
    setDefault("coarsestTolerance", coupledDefaults::coarsestTolerance);
    setDefault("coarsestMaxIter", coupledDefaults::coarsestMaxIter);
    setDefault("cacheAgglomeration", coupledDefaults::cacheAgglomeration);
    dict_.set("cycleType", cycleName(cycle_));

    agglomerate();

    const GAMGAgglomeration& agg = *aggPtr_;
    const label nCoarse = agg.size();

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

    const label nLev = nLevels();
    smoothers_.resize(nLev - 1);
    for (label l = 0; l < nLev - 1; ++l)
    {
        smoothers_.set(l, blockSmoother::New(matrixLevel(l), dict_).ptr());
    }

    // Coarsest level (6.3.3): dense LU if small and on one rank, else Krylov
    {
        const blockLduMatrix4& Ac = matrixLevel(L());
        const bool oneRank =
            !UPstream::parRun() || UPstream::nProcs(Ac.comm()) == 1;
        useDenseLU_ =
            oneRank && Ac.interfaces().empty()
         && Ac.nCells() <= denseLUMaxCells_;

        if (!useDenseLU_)
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
            coarsestSolver_ = blockSolver::New(Ac, cd);
        }
        Info<< "blockGAMG: coarsest level " << L() << " ("
            << returnReduce(Ac.nCells(), sumOp<label>()) << " cells): "
            << (useDenseLU_ ? "dense LU" : dict_.get<word>("coarsestSolver"))
            << endl;
    }

    allocateWork();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockGAMG::measure
(
    const GAMGAgglomeration& agg,
    reduceScalar& cop,
    List<reduceScalar>& ratios
) const
{
    const label nC = agg.size();

    // Per level: cells and nnz blocks, summed over ranks in one call
    List<reduceScalar> cells(nC + 1, Zero);
    List<reduceScalar> nnz(nC + 1, Zero);
    cells[0] = reduceScalar(fine_.nCells());
    nnz[0] = reduceScalar(fine_.nnzBlocks());

    for (label l = 0; l < nC; ++l)
    {
        cells[l + 1] = reduceScalar(agg.nCells(l));
        reduceScalar n = cells[l + 1] + 2*reduceScalar(agg.nFaces(l));

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
        nnz[l + 1] = n;
    }

    doubleReduce::parSum(cells.data(), cells.size(), fine_.comm());
    doubleReduce::parSum(nnz.data(), nnz.size(), fine_.comm());

    reduceScalar total = 0;
    for (const reduceScalar v : nnz)
    {
        total += v;
    }
    // GUARD: a mesh has at least one cell globally
    cop = total/std::max(nnz[0], reduceScalar(1));

    ratios.resize(nC);
    for (label l = 0; l < nC; ++l)
    {
        // GUARD: coarse levels are never empty globally
        ratios[l] = cells[l]/std::max(cells[l + 1], reduceScalar(1));
    }
}


bool Foam::blockGAMG::ratiosOk(const List<reduceScalar>& ratios) const
{
    for (const reduceScalar r : ratios)
    {
        if (r < minRatio_)
        {
            return false;
        }
    }
    return true;
}


void Foam::blockGAMG::agglomerate()
{
    const lduMesh& mesh = fine_.mesh();
    const label m0 = dict_.get<label>("mergeLevels");
    const label mMax = max(m0, maxMergeLevels_);

    reduceScalar cop = 0;
    List<reduceScalar> ratios;

    for (label m = m0; m <= mMax; ++m)
    {
        dictionary aggDict(dict_);
        aggDict.set("mergeLevels", m);
        aggDict.set("name", word("blockGAMGAgglomeration_m" + Foam::name(m)));

        const GAMGAgglomeration& agg = GAMGAgglomeration::New(mesh, aggDict);
        measure(agg, cop, ratios);

        const bool copOk = cop <= maxCop_;
        const bool ratioOk = !needsRatio(cycle_) || ratiosOk(ratios);

        Info<< "blockGAMG: mergeLevels " << m << " levels " << agg.size() + 1
            << " C_op " << cop << " ratios " << flatOutput(ratios)
            << (copOk && ratioOk ? " accepted" : " rejected") << endl;

        if (copOk && ratioOk)
        {
            aggPtr_ = &agg;
            mergeLevelsUsed_ = m;
            Cop_ = cop;
            ratios_ = ratios;
            return;
        }
    }

    FatalErrorInFunction
        << "No agglomeration with mergeLevels " << m0 << ".." << mMax
        << " satisfies maxOperatorComplexity " << maxCop_
        << (needsRatio(cycle_)
            ? std::string(" and the coarsening-ratio rule r_l >= ")
              + Foam::name(minRatio_) + " (cycleType " + cycleName(cycle_)
              + ", 6.3.1)"
            : std::string())
        << nl << "Last measured: C_op " << cop << ", ratios "
        << flatOutput(ratios) << exit(FatalError);
}


void Foam::blockGAMG::allocateWork()
{
    const label nLev = nLevels();
    for (auto* v : {&b_, &e_, &r_, &t_, &e2_, &z1_, &q1_, &z2_, &q2_, &rk_})
    {
        v->resize(nLev);
    }
    for (label l = 0; l < nLev; ++l)
    {
        const label n = matrixLevel(l).nRows();
        r_[l].resize(n, Zero);
        if (l > 0)
        {
            b_[l].resize(n, Zero);
            e_[l].resize(n, Zero);
            t_[l].resize(n, Zero);
            e2_[l].resize(n, Zero);
            if (cycle_ == cycleKind::K)
            {
                z1_[l].resize(n, Zero);
                q1_[l].resize(n, Zero);
                z2_[l].resize(n, Zero);
                q2_[l].resize(n, Zero);
                rk_[l].resize(n, Zero);
            }
        }
    }
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


void Foam::blockGAMG::restrictVector
(
    const label l,
    const blockScalarUList& r,
    blockScalarUList& bc
) const
{
    const labelField& restrictAddr = aggPtr_->restrictAddressing(l);
    bc = Zero;
    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockDim; ++k)
        {
            bc[c*blockDim + k] += r[celli*blockDim + k];
        }
    }
}


void Foam::blockGAMG::prolongAdd
(
    const label l,
    const blockScalarUList& e,
    blockScalarUList& x
) const
{
    const labelField& restrictAddr = aggPtr_->restrictAddressing(l);
    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockDim; ++k)
        {
            x[celli*blockDim + k] += e[c*blockDim + k];
        }
    }
}


void Foam::blockGAMG::factoriseDense()
{
    if (!useDenseLU_)
    {
        return;
    }

    const blockLduMatrix4& A = matrixLevel(L());
    const label nc = A.nCells();
    const label n = blockDim*nc;

    denseLU_.resize_nocopy(n*n);
    denseLU_ = Zero;
    densePivot_.resize_nocopy(n);

    auto addBlock = [&](const label rc, const label cc, const blockScalar* B)
    {
        for (label r = 0; r < blockDim; ++r)
        {
            for (label c = 0; c < blockDim; ++c)
            {
                denseLU_[(rc*blockDim + r)*n + cc*blockDim + c] +=
                    toDouble(B[r*blockDim + c]);
            }
        }
    };

    for (label celli = 0; celli < nc; ++celli)
    {
        addBlock(celli, celli, A.diagBlock(celli));
    }
    const labelUList& lAddr = A.lduAddr().lowerAddr();
    const labelUList& uAddr = A.lduAddr().upperAddr();
    forAll(lAddr, facei)
    {
        addBlock(lAddr[facei], uAddr[facei], A.upper().cdata() + facei*blockSize);
        addBlock(uAddr[facei], lAddr[facei], A.lower().cdata() + facei*blockSize);
    }

    // LU with partial pivoting (in place), pivot guard as block4Ops
    const reduceScalar guard =
        dict_.getOrDefault<doubleScalar>("pivotGuard", coupledDefaults::pivotGuard);

    for (label k = 0; k < n; ++k)
    {
        label piv = k;
        reduceScalar big = std::abs(denseLU_[k*n + k]);
        for (label i = k + 1; i < n; ++i)
        {
            const reduceScalar a = std::abs(denseLU_[i*n + k]);
            if (a > big)
            {
                big = a;
                piv = i;
            }
        }
        densePivot_[k] = piv;
        if (piv != k)
        {
            for (label j = 0; j < n; ++j)
            {
                std::swap(denseLU_[k*n + j], denseLU_[piv*n + j]);
            }
        }
        // GUARD: pivot guard
        if (std::abs(denseLU_[k*n + k]) < guard)
        {
            denseLU_[k*n + k] += (denseLU_[k*n + k] < 0 ? -guard : guard);
        }
        const reduceScalar rp = 1.0/denseLU_[k*n + k];  // GUARD: guarded
        for (label i = k + 1; i < n; ++i)
        {
            const reduceScalar f = (denseLU_[i*n + k] *= rp);
            if (f != 0)
            {
                for (label j = k + 1; j < n; ++j)
                {
                    denseLU_[i*n + j] -= f*denseLU_[k*n + j];
                }
            }
        }
    }
}


void Foam::blockGAMG::solveCoarsest
(
    blockScalarUList& x,
    const blockScalarUList& b
) const
{
    if (!useDenseLU_)
    {
        x = Zero;
        const blockSolverPerformance perf = coarsestSolver_->solve(x, b);
        nCoarsestIters_ = perf.nIterations;
        return;
    }

    const label n = b.size();
    reduceScalarList y(n);
    for (label i = 0; i < n; ++i)
    {
        y[i] = toDouble(b[i]);
    }
    for (label k = 0; k < n; ++k)
    {
        const label p = densePivot_[k];
        if (p != k)
        {
            std::swap(y[k], y[p]);
        }
    }
    for (label i = 0; i < n; ++i)
    {
        reduceScalar s = y[i];
        for (label j = 0; j < i; ++j)
        {
            s -= denseLU_[i*n + j]*y[j];
        }
        y[i] = s;
    }
    for (label i = n - 1; i >= 0; --i)
    {
        reduceScalar s = y[i];
        for (label j = i + 1; j < n; ++j)
        {
            s -= denseLU_[i*n + j]*y[j];
        }
        y[i] = s/denseLU_[i*n + i];  // GUARD: pivots guarded
    }
    for (label i = 0; i < n; ++i)
    {
        x[i] = narrow(y[i]);
    }
    nCoarsestIters_ = 0;
}


void Foam::blockGAMG::cycle
(
    const label l,
    const blockScalarUList& b,
    blockScalarUList& x,
    const cycleKind kind
) const
{
    if (l == L())
    {
        solveCoarsest(x, b);
        return;
    }

    const blockLduMatrix4& A = matrixLevel(l);
    const label nPre = (l == 0 ? nFinestSweeps_ : nPreSweeps_);
    const label nPost = (l == 0 ? nFinestSweeps_ : nPostSweeps_);

    x = Zero;
    if (nPre > 0)
    {
        smoothers_[l].smooth(x, b, nPre);
    }

    A.residual(r_[l], x, b);
    restrictVector(l, r_[l], b_[l + 1]);

    blockScalarList& ec = e_[l + 1];

    switch (kind)
    {
        case cycleKind::V:
        {
            cycle(l + 1, b_[l + 1], ec, cycleKind::V);
            break;
        }
        case cycleKind::F:
        {
            cycle(l + 1, b_[l + 1], ec, cycleKind::F);
            if (l + 1 < L())
            {
                matrixLevel(l + 1).residual(t_[l + 1], ec, b_[l + 1]);
                cycle(l + 1, t_[l + 1], e2_[l + 1], cycleKind::V);
                const blockScalarList& e2 = e2_[l + 1];
                forAll(ec, i)
                {
                    ec[i] += e2[i];
                }
            }
            break;
        }
        case cycleKind::W:
        {
            cycle(l + 1, b_[l + 1], ec, cycleKind::W);
            if (l + 1 < L())
            {
                matrixLevel(l + 1).residual(t_[l + 1], ec, b_[l + 1]);
                cycle(l + 1, t_[l + 1], e2_[l + 1], cycleKind::W);
                const blockScalarList& e2 = e2_[l + 1];
                forAll(ec, i)
                {
                    ec[i] += e2[i];
                }
            }
            break;
        }
        case cycleKind::K:
        {
            kstep(l + 1, b_[l + 1], ec);
            break;
        }
    }

    prolongAdd(l, ec, x);

    if (nPost > 0)
    {
        smoothers_[l].smooth(x, b, nPost);
    }
}


void Foam::blockGAMG::kstep
(
    const label l,
    const blockScalarUList& b,
    blockScalarUList& e
) const
{
    if (l == L())
    {
        solveCoarsest(e, b);
        return;
    }

    const blockLduMatrix4& A = matrixLevel(l);
    const label comm = A.comm();
    const label n = A.nRows();

    blockScalarList& z1 = z1_[l];
    blockScalarList& q1 = q1_[l];
    blockScalarList& z2 = z2_[l];
    blockScalarList& q2 = q2_[l];
    blockScalarList& r1 = rk_[l];

    const reduceScalar r0norm = doubleReduce::norm2(b, comm);

    // First GCR step: z1 = M^-1 b, q1 = A z1
    cycle(l, b, z1, cycleKind::K);
    A.Amul(q1, z1);

    const FixedList<reduceScalar, 2> d1 = doubleReduce::dot2(q1, b, q1, q1, comm);
    // GUARD: <q1,q1> > 0 unless the correction vanishes
    const reduceScalar a1 = d1[0]/std::max(d1[1], doubleScalarVSMALL);
    const blockScalar a1f = narrow(a1);

    for (label i = 0; i < n; ++i)
    {
        e[i] = a1f*z1[i];
        r1[i] = b[i] - a1f*q1[i];
    }

    if (kMaxSteps_ < 2)
    {
        return;
    }
    const reduceScalar r1norm = doubleReduce::norm2(r1, comm);
    if (r1norm <= kThreshold_*r0norm)
    {
        return;
    }

    // Second GCR step, orthogonalised against q1
    cycle(l, r1, z2, cycleKind::K);
    A.Amul(q2, z2);

    // coefficient from the pre-orthogonalisation q2 (6.3.2)
    const reduceScalar c =
        doubleReduce::dot(q2, q1, comm)/std::max(d1[1], doubleScalarVSMALL);
    const blockScalar cf = narrow(c);
    for (label i = 0; i < n; ++i)
    {
        q2[i] -= cf*q1[i];
        z2[i] -= cf*z1[i];
    }

    const FixedList<reduceScalar, 2> d2 = doubleReduce::dot2(q2, r1, q2, q2, comm);
    // GUARD: <q2,q2> > 0 unless the correction vanishes
    const reduceScalar a2 = d2[0]/std::max(d2[1], doubleScalarVSMALL);
    const blockScalar a2f = narrow(a2);
    for (label i = 0; i < n; ++i)
    {
        e[i] += a2f*z2[i];
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


bool Foam::blockGAMG::setCycleType(const cycleKind c) const
{
    if (needsRatio(c) && !ratiosOk(ratios_))
    {
        return false;
    }
    if (c == cycleKind::K && z1_.size() > 1 && z1_[1].empty())
    {
        // K workspace was not allocated for the start-up cycle type
        for (label l = 1; l < nLevels(); ++l)
        {
            const label n = matrixLevel(l).nRows();
            z1_[l].resize(n, Zero);
            q1_[l].resize(n, Zero);
            z2_[l].resize(n, Zero);
            q2_[l].resize(n, Zero);
            rk_[l].resize(n, Zero);
        }
    }
    cycle_ = c;
    return true;
}


void Foam::blockGAMG::update()
{
    if (!dict_.get<bool>("cacheAgglomeration"))
    {
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

    factoriseDense();
}


void Foam::blockGAMG::apply
(
    blockScalarUList& x,
    const blockScalarUList& b
) const
{
    cycle(0, b, x, cycle_);
}


void Foam::blockGAMG::writeStats(Ostream& os) const
{
    const labelList n = globalCellsPerLevel();
    os  << "blockGAMG: levels " << nLevels()
        << ", mergeLevels " << mergeLevelsUsed_
        << ", C_op " << Cop_
        << ", cycle " << cycleName(cycle_) << nl
        << "blockGAMG: cells per level " << flatOutput(n) << nl
        << "blockGAMG: coarsening ratios " << flatOutput(ratios_) << endl;
}


Foam::dictionary Foam::blockGAMG::statsDict() const
{
    dictionary d;
    d.add("nLevels", nLevels());
    d.add("mergeLevels", mergeLevelsUsed_);
    d.add("Cop", Cop_);
    d.add("cycleType", cycleName(cycle_));
    d.add("nPostSweeps", nPostSweeps_);
    d.add("cellsPerLevel", globalCellsPerLevel());
    d.add("ratios", List<doubleScalar>(ratios_));
    d.add("denseCoarsest", useDenseLU_);
    return d;
}


void Foam::blockGAMG::writeSettings(dictionary& dict) const
{
    dict.merge(dict_);
    dict.set("cycleType", cycleName(cycle_));
    dict.set("kCycleThreshold", kThreshold_);
    dict.set("kCycleMaxSteps", kMaxSteps_);
    dict.set("nPreSweeps", nPreSweeps_);
    dict.set("nPostSweeps", nPostSweeps_);
    dict.set("nFinestSweeps", nFinestSweeps_);
    dict.set("maxOperatorComplexity", maxCop_);
    dict.set("minCoarseningRatio", minRatio_);
    dict.set("maxMergeLevels", maxMergeLevels_);
    dict.set("mergeLevelsUsed", mergeLevelsUsed_);
    dict.set("denseLUMaxCells", denseLUMaxCells_);
}


// ************************************************************************* //
