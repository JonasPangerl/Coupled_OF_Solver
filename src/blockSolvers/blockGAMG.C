/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGAMG.H"
#include "blockGAMGProcAgglomeration.H"
#include "blockPreconditioner.H"
#include "block4Ops.H"
#include "doubleReduce.H"
#include "coupledDefaults.H"
#include "lduMesh.H"
#include "lduPrimitiveMesh.H"
#include "globalIndex.H"
#include "PstreamBuffers.H"
#include "objectRegistry.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "pairGAMGAgglomeration.H"
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdint>

// * * * * * * * * * * * * blockGAMGProcAgglomeration  * * * * * * * * * * * //
// Implemented here (not in its own translation unit) so that the library
// file list does not change.

namespace Foam
{
    defineTypeNameAndDebug(blockGAMGProcAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGProcAgglomeration,
        blockGAMGProcAgglomeration,
        GAMGAgglomeration
    );
}


Foam::blockGAMGProcAgglomeration::blockGAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    GAMGProcAgglomeration(agglom, controlDict),
    cellsPerRank_
    (
        controlDict.getOrDefault<label>
        (
            "procAgglomCellsPerRank",
            coupledDefaults::procAgglomCellsPerRank
        )
    ),
    divisor_
    (
        controlDict.getOrDefault<label>
        (
            "procAgglomDivisor",
            coupledDefaults::procAgglomDivisor
        )
    )
{
    if (cellsPerRank_ < 1 || divisor_ < 1)
    {
        FatalIOErrorInFunction(controlDict)
            << "procAgglomCellsPerRank (" << cellsPerRank_
            << ") and procAgglomDivisor (" << divisor_
            << ") must be >= 1" << exit(FatalIOError);
    }
}


Foam::label Foam::blockGAMGProcAgglomeration::targetRanks
(
    const reduceScalar nCells,
    const label nRanks,
    const label cellsPerRank,
    const label divisor
)
{
    if (nCells < reduceScalar(cellsPerRank))
    {
        return 1;
    }
    if (nCells < reduceScalar(cellsPerRank)*reduceScalar(nRanks))
    {
        return max(label(1), nRanks/divisor);
    }
    return nRanks;
}


Foam::labelList Foam::blockGAMGProcAgglomeration::groupMap
(
    const label nProcs,
    const label nGroups
)
{
    labelList map(nProcs);
    forAll(map, proci)
    {
        // 64-bit product: no overflow for any rank count
        map[proci] = label
        (
            (int64_t(proci)*int64_t(nGroups))/int64_t(nProcs)
        );
    }
    return map;
}


bool Foam::blockGAMGProcAgglomeration::agglomerate()
{
    const label nCoarse = agglom_.size();
    const label comm0 = agglom_.mesh().comm();
    const label nRanks = UPstream::nProcs(comm0);

    // Levels 1..nCoarse-1 can be gathered natively (level 0 never, the
    // coarsest level nCoarse is gathered by blockGAMG::gatherCoarsest)
    if (nCoarse < 2 || nRanks < 2)
    {
        return true;
    }

    // Global cell count per level. At this point (end of the local
    // agglomeration, before any processor agglomeration) every rank holds
    // every level.
    reduceScalarList cells(nCoarse + 1, Zero);
    for (label l = 1; l <= nCoarse; ++l)
    {
        cells[l] = reduceScalar(agglom_.nCells(l - 1));
    }
    doubleReduce::parSum(cells.data(), cells.size(), comm0);

    label current = nRanks;

    for (label l = 1; l < nCoarse; ++l)
    {
        const label target =
            targetRanks(cells[l], nRanks, cellsPerRank_, divisor_);

        if (target >= current)
        {
            continue;
        }

        // Ranks agglomerated away at an earlier step do not hold the level
        // and take no part (as native manualGAMGProcAgglomeration)
        if (agglom_.hasMeshLevel(l))
        {
            const label levelComm = agglom_.meshLevel(l).comm();
            const label nProcs = UPstream::nProcs(levelComm);

            if (nProcs > 1 && UPstream::myProcNo(levelComm) != -1)
            {
                const labelList procAgglomMap(groupMap(nProcs, target));

                labelList masterProcs;
                List<label> agglomProcIDs;
                GAMGAgglomeration::calculateRegionMaster
                (
                    levelComm,
                    procAgglomMap,
                    masterProcs,
                    agglomProcIDs
                );

                // Communicator for the processor-agglomerated level
                comms_.push_back
                (
                    UPstream::newCommunicator(levelComm, masterProcs)
                );

                // Native gathering of the level and re-agglomeration of the
                // coarser levels on the masters
                GAMGProcAgglomeration::agglomerate
                (
                    l,
                    procAgglomMap,
                    masterProcs,
                    agglomProcIDs,
                    comms_.back()
                );
            }
        }

        current = target;
    }

    return true;
}


// * * * * * * * * * * * * blockPairAgglomeration  * * * * * * * * * * * * * //

namespace Foam
{

//- Native pair agglomeration on given finest-level face weights
//  (precond-research: weights from the block matrix, see
//  blockGAMG::matrixFaceWeights)
class blockPairAgglomeration
:
    public pairGAMGAgglomeration
{
public:

    blockPairAgglomeration
    (
        const lduMesh& mesh,
        const dictionary& controlDict,
        const tmp<scalarField>& tweights
    )
    :
        pairGAMGAgglomeration(mesh, controlDict)
    {
        agglomerate(nCellsInCoarsestLevel_, 0, tweights(), true);
    }
};

} // End namespace Foam


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
    scaleMode_(0),
    aggWeights_
    (
        dict.getOrDefault<word>
        (
            "agglomerationWeights",
            word(coupledDefaults::agglomerationWeights)
        )
    ),
    reaggInterval_
    (
        dict.getOrDefault<label>
        (
            "reagglomerateInterval",
            coupledDefaults::reagglomerateInterval
        )
    ),
    nUpdates_(0),
    aggFromMatrix_(false),
    procAgglomType_
    (
        dict.getOrDefault<word>("processorAgglomerator", "masterCoarsest")
    ),
    procAgglomCellsPerRank_
    (
        dict.getOrDefault<label>
        (
            "procAgglomCellsPerRank",
            coupledDefaults::procAgglomCellsPerRank
        )
    ),
    procAgglomDivisor_
    (
        dict.getOrDefault<label>
        (
            "procAgglomDivisor",
            coupledDefaults::procAgglomDivisor
        )
    ),
    aggPtr_(nullptr),
    mergeLevelsUsed_(-1),
    Cop_(0),
    ratios_(),
    L_(0),
    cellsPerLevel_(),
    ranksPerLevel_(),
    levelMesh_(),
    gather_(),
    coarsestMesh_(),
    coarsestAllComm_(-1),
    coarsestAgglomComm_(-1),
    coarse_(),
    smoothers_(),
    coarsestSolver_(),
    useDenseLU_(false),
    denseLU_(),
    densePivot_(),
    nCoarsestIters_(0),
    hierarchyVersion_(0),
    diag_(nullptr)
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
    {
        const word sc
        (
            dict.getOrDefault<word>
            (
                "scaleCorrection",
                word(coupledDefaults::scaleCorrection)
            )
        );
        if (sc == "none") scaleMode_ = 0;
        else if (sc == "finest") scaleMode_ = 1;
        else if (sc == "all") scaleMode_ = 2;
        else
        {
            FatalIOErrorInFunction(dict)
                << "scaleCorrection " << sc << ": valid none finest all"
                << exit(FatalIOError);
        }
        dict_.set("scaleCorrection", sc);
    }
    if
    (
        aggWeights_ != "geometric" && aggWeights_ != "momentum"
     && aggWeights_ != "pressure" && aggWeights_ != "combined"
    )
    {
        FatalIOErrorInFunction(dict)
            << "agglomerationWeights " << aggWeights_
            << ": valid geometric momentum pressure combined"
            << exit(FatalIOError);
    }
    dict_.set("agglomerationWeights", aggWeights_);
    dict_.set("reagglomerateInterval", reaggInterval_);
    dict_.set("cycleType", cycleName(cycle_));
    dict_.set("processorAgglomerator", procAgglomType_);
    if (ruleMode())
    {
        if (procAgglomCellsPerRank_ < 1 || procAgglomDivisor_ < 1)
        {
            FatalIOErrorInFunction(dict)
                << "procAgglomCellsPerRank (" << procAgglomCellsPerRank_
                << ") and procAgglomDivisor (" << procAgglomDivisor_
                << ") must be >= 1" << exit(FatalIOError);
        }
        dict_.set("procAgglomCellsPerRank", procAgglomCellsPerRank_);
        dict_.set("procAgglomDivisor", procAgglomDivisor_);
    }

    buildHierarchy();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockGAMG::~blockGAMG()
{
    clearHierarchy();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blockGAMG::buildHierarchy()
{
    ++hierarchyVersion_;
    agglomerate();

    const GAMGAgglomeration& agg = *aggPtr_;
    L_ = agg.size();

    // Levels held by this rank
    levelMesh_.clear();
    levelMesh_.resize(L_ + 1);
    levelMesh_.set(0, &fine_.mesh());
    for (label l = 1; l <= L_; ++l)
    {
        if (agg.hasMeshLevel(l))
        {
            levelMesh_.set(l, &agg.meshLevel(l));
        }
    }

    // Processor agglomeration (6.3.4)
    gather_.clear();
    gather_.resize(L_ + 1);
    setNativeGathering();
    gatherCoarsest();

    // Vector offsets of the gathered levels: blockDim values per cell
    forAll(gather_, l)
    {
        if (gather_.set(l))
        {
            gatherLevel& g = gather_[l];
            g.vecOffsets.resize_nocopy(g.cellOffsets.size());
            forAll(g.vecOffsets, i)
            {
                g.vecOffsets[i] = blockDim*g.cellOffsets[i];
            }
        }
    }

    // Ranks per level (global)
    const label nLevG = cellsPerLevel_.size();
    ranksPerLevel_.resize_nocopy(nLevG);
    for (label l = 0; l < nLevG; ++l)
    {
        ranksPerLevel_[l] = (hasLevel(l) ? 1 : 0);
    }
    Foam::reduce
    (
        ranksPerLevel_.data(),
        int(ranksPerLevel_.size()),
        sumOp<label>(),
        UPstream::msgType(),
        fine_.comm()
    );

    // Coarse matrices on the levels held
    coarse_.clear();
    coarse_.resize(L_);
    for (label l = 1; l <= L_; ++l)
    {
        if (hasLevel(l))
        {
            coarse_.set(l - 1, new blockLduMatrix4(levelMesh_[l]));
        }
    }

    smoothers_.clear();
    smoothers_.resize(L_);
    for (label l = 0; l < L_; ++l)
    {
        if (hasLevel(l))
        {
            smoothers_.set(l, blockSmoother::New(matrixLevel(l), dict_).ptr());
        }
    }

    // Coarsest level (6.3.3): dense LU if small and on one rank, else Krylov
    useDenseLU_ = false;
    coarsestSolver_.reset(nullptr);
    if (hasLevel(L_))
    {
        const blockLduMatrix4& Ac = matrixLevel(L_);
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
    }
    {
        const bool anyDense =
            returnReduceOr(useDenseLU_, fine_.comm());
        Info<< "blockGAMG: coarsest level " << nLevG - 1 << " ("
            << cellsPerLevel_.last() << " cells, "
            << ranksPerLevel_.last() << " rank(s)): "
            << (anyDense ? "dense LU" : dict_.get<word>("coarsestSolver"))
            << endl;
    }

    allocateWork();
}


void Foam::blockGAMG::clearHierarchy()
{
    coarsestSolver_.reset(nullptr);
    smoothers_.clear();
    coarse_.clear();
    useDenseLU_ = false;
    denseLU_.clear();
    densePivot_.clear();
    gather_.clear();
    levelMesh_.clear();
    coarsestMesh_.reset(nullptr);

    // Children before parents
    if (coarsestAgglomComm_ >= 0)
    {
        UPstream::freeCommunicator(coarsestAgglomComm_);
        coarsestAgglomComm_ = -1;
    }
    if (coarsestAllComm_ >= 0)
    {
        UPstream::freeCommunicator(coarsestAllComm_);
        coarsestAllComm_ = -1;
    }
}


void Foam::blockGAMG::measure
(
    const GAMGAgglomeration& agg,
    reduceScalar& cop,
    List<reduceScalar>& ratios,
    labelList& cellsOut
) const
{
    const label comm = fine_.comm();
    const label nC = agg.size();

    // Number of levels: the same on all ranks except for native
    // agglomerators that add master-only levels
    const label nLevG =
        returnReduce(nC, maxOp<label>(), UPstream::msgType(), comm) + 1;

    // Per level: cells and nnz blocks of the levels held by this rank,
    // summed over ranks in one call. The sums do not depend on processor
    // agglomeration (a merged processor face counts twice before and after).
    List<reduceScalar> cells(nLevG, Zero);
    List<reduceScalar> nnz(nLevG, Zero);
    cells[0] = reduceScalar(fine_.nCells());
    nnz[0] = reduceScalar(fine_.nnzBlocks());

    for (label l = 1; l <= nC; ++l)
    {
        if (!agg.hasMeshLevel(l))
        {
            continue;
        }
        const lduMesh& m = agg.meshLevel(l);
        cells[l] = reduceScalar(m.lduAddr().size());
        reduceScalar n =
            cells[l] + 2*reduceScalar(m.lduAddr().lowerAddr().size());

        const lduInterfacePtrsList ifaces = m.interfaces();
        forAll(ifaces, i)
        {
            if
            (
                ifaces.set(i)
             && blockLduInterface::isBlockCoupled(ifaces[i])
            )
            {
                n += reduceScalar(ifaces[i].faceCells().size());
            }
        }
        nnz[l] = n;
    }

    doubleReduce::parSum(cells.data(), cells.size(), comm);
    doubleReduce::parSum(nnz.data(), nnz.size(), comm);

    reduceScalar total = 0;
    for (const reduceScalar v : nnz)
    {
        total += v;
    }
    // GUARD: a mesh has at least one cell globally
    cop = total/std::max(nnz[0], reduceScalar(1));

    ratios.resize(nLevG - 1);
    for (label l = 0; l < nLevG - 1; ++l)
    {
        // GUARD: coarse levels are never empty globally
        ratios[l] = cells[l]/std::max(cells[l + 1], reduceScalar(1));
    }

    cellsOut.resize_nocopy(nLevG);
    forAll(cellsOut, l)
    {
        cellsOut[l] = label(std::llround(cells[l]));
    }
}


Foam::tmp<Foam::scalarField> Foam::blockGAMG::matrixFaceWeights() const
{
    // Strength of the coupling across each internal face, from the
    // finest block matrix:
    //   momentum  max(|U_f(0,0)|, |L_f(0,0)|)            (u-u block)
    //   pressure  max(|U_f(3,3)|, |L_f(3,3)|)            (p-p block)
    //   combined  sum over (0,0) and (3,3) of the classical strength
    //             max(|U|,|L|)/sqrt(|D_P| |D_N|)
    const label uu = 0;
    const label pp = blockSize - 1;
    const lduAddressing& addr = fine_.lduAddr();
    const labelUList& lAddr = addr.lowerAddr();
    const labelUList& uAddr = addr.upperAddr();
    const blockScalarList& fU = fine_.upper();
    const blockScalarList& fL = fine_.lower();
    const blockScalarList& fD = fine_.diag();

    auto strength = [&](const label facei, const label k)
    {
        return std::max
        (
            std::abs(toDouble(fU[facei*blockSize + k])),
            std::abs(toDouble(fL[facei*blockSize + k]))
        );
    };
    auto normalised = [&](const label facei, const label k)
    {
        const reduceScalar dP =
            std::abs(toDouble(fD[lAddr[facei]*blockSize + k]));
        const reduceScalar dN =
            std::abs(toDouble(fD[uAddr[facei]*blockSize + k]));
        // GUARD: a zero diagonal gives a zero-strength face
        return strength(facei, k)
            /std::max(std::sqrt(dP*dN), doubleScalarVSMALL);
    };

    auto tw = tmp<scalarField>::New(lAddr.size(), Zero);
    scalarField& w = tw.ref();
    forAll(w, facei)
    {
        if (aggWeights_ == "momentum")
        {
            w[facei] = strength(facei, uu);
        }
        else if (aggWeights_ == "pressure")
        {
            w[facei] = strength(facei, pp);
        }
        else
        {
            w[facei] = normalised(facei, uu) + normalised(facei, pp);
        }
    }
    return tw;
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
    labelList cells;

    for (label m = m0; m <= mMax; ++m)
    {
        dictionary aggDict(dict_);
        aggDict.set("mergeLevels", m);

        // Processor agglomeration (6.3.4)
        if (procAgglomType_ == "none")
        {
            aggDict.remove("processorAgglomerator");
        }
        else if (ruleMode())
        {
            aggDict.set
            (
                "processorAgglomerator",
                word(blockGAMGProcAgglomeration::typeName)
            );
        }
        else if (procAgglomType_ == "nativeMasterCoarsest")
        {
            aggDict.set("processorAgglomerator", word("masterCoarsest"));
        }
        // else: native type passed verbatim

        aggDict.set
        (
            "name",
            word
            (
                "blockGAMGAgglomeration_m" + Foam::name(m)
              + "_" + procAgglomType_
            )
        );

        const GAMGAgglomeration* aggP = nullptr;
        if (aggFromMatrix_)
        {
            // Pair agglomeration on block-matrix face weights; replaces a
            // previous one of the same name (re-agglomeration)
            const word nm
            (
                aggDict.get<word>("name") + "_" + aggWeights_
            );
            aggDict.set("name", nm);
            const GAMGAgglomeration* old =
                mesh.thisDb().cfindObject<GAMGAgglomeration>(nm);
            if (old)
            {
                mesh.thisDb().checkOut(const_cast<GAMGAgglomeration*>(old));
            }
            autoPtr<GAMGAgglomeration> p
            (
                new blockPairAgglomeration(mesh, aggDict, matrixFaceWeights())
            );
            aggP = &regIOobject::store(p);
        }
        else
        {
            aggP = &GAMGAgglomeration::New(mesh, aggDict);
        }
        const GAMGAgglomeration& agg = *aggP;
        measure(agg, cop, ratios, cells);

        const bool copOk = cop <= maxCop_;
        const bool ratioOk = !needsRatio(cycle_) || ratiosOk(ratios);

        Info<< "blockGAMG: mergeLevels " << m << " levels " << cells.size()
            << " C_op " << cop << " ratios " << flatOutput(ratios)
            << (copOk && ratioOk ? " accepted" : " rejected") << endl;

        if (copOk && ratioOk)
        {
            aggPtr_ = &agg;
            mergeLevelsUsed_ = m;
            Cop_ = cop;
            ratios_ = ratios;
            cellsPerLevel_ = cells;
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


void Foam::blockGAMG::setNativeGathering()
{
    const GAMGAgglomeration& agg = *aggPtr_;

    if (!agg.processorAgglomerate())
    {
        return;
    }

    for (label l = 1; l <= L_; ++l)
    {
        // Level l is gathered from the local coarse operators of level l-1
        if (!hasLevel(l - 1) || !agg.hasProcMesh(l))
        {
            continue;
        }

        auto* gPtr = new gatherLevel;
        gatherLevel& g = *gPtr;
        g.comm = agg.agglomCommunicator(l);
        g.master = (UPstream::myProcNo(g.comm) == 0);
        g.procIDs = identity(UPstream::nProcs(g.comm));

        if (g.master)
        {
            g.cellOffsets = agg.cellOffsets(l);
            g.faceMap = agg.faceMap(l);
            g.boundaryMap = agg.boundaryMap(l);
            g.boundaryFaceMap = agg.boundaryFaceMap(l);
        }
        else if (hasLevel(l))
        {
            FatalErrorInFunction
                << "Rank holds processor-agglomerated level " << l
                << " without being the master of its group"
                << abort(FatalError);
        }

        gather_.set(l, gPtr);
    }
}


void Foam::blockGAMG::gatherCoarsest()
{
    if (!ruleMode() || L_ < 1 || !UPstream::parRun())
    {
        return;
    }

    const label nRanks = UPstream::nProcs(fine_.comm());
    const label target = blockGAMGProcAgglomeration::targetRanks
    (
        reduceScalar(cellsPerLevel_.last()),
        nRanks,
        procAgglomCellsPerRank_,
        procAgglomDivisor_
    );

    const label current = returnReduce
    (
        label(hasLevel(L_) ? 1 : 0),
        sumOp<label>(),
        UPstream::msgType(),
        fine_.comm()
    );

    if (target >= current || !hasLevel(L_))
    {
        return;
    }

    // This rank holds the coarsest level, which is to be gathered onto
    // 'target' masters. Same native steps as the native processor
    // agglomeration of a level (masterCoarsestGAMGProcAgglomeration +
    // GAMGAgglomeration::procAgglomerateLduAddressing), but without the
    // restriction addressing below the level, which does not exist.

    const lduMesh& localMesh = levelMesh_[L_];
    const label levelComm = localMesh.comm();
    const label nProcs = UPstream::nProcs(levelComm);

    if (nProcs != current)
    {
        FatalErrorInFunction
            << "Coarsest level communicator has " << nProcs
            << " ranks, expected " << current << abort(FatalError);
    }

    const labelList procAgglomMap
    (
        blockGAMGProcAgglomeration::groupMap(nProcs, target)
    );

    labelList masterProcs;
    List<label> agglomProcIDs;
    GAMGAgglomeration::calculateRegionMaster
    (
        levelComm,
        procAgglomMap,
        masterProcs,
        agglomProcIDs
    );

    // Communicator of the gathered level (the masters) and of my group
    coarsestAllComm_ = UPstream::newCommunicator(levelComm, masterProcs);
    coarsestAgglomComm_ = UPstream::newCommunicator(levelComm, agglomProcIDs);

    // Collect the meshes of the group on its master
    PtrList<lduPrimitiveMesh> otherMeshes;
    lduPrimitiveMesh::gather(coarsestAgglomComm_, localMesh, otherMeshes);

    auto* gPtr = new gatherLevel;
    gatherLevel& g = *gPtr;
    g.comm = coarsestAgglomComm_;
    g.master = (UPstream::myProcNo(g.comm) == 0);
    g.procIDs = identity(UPstream::nProcs(g.comm));

    if (g.master)
    {
        labelList faceOffsets;
        coarsestMesh_.reset
        (
            new lduPrimitiveMesh
            (
                coarsestAllComm_,
                procAgglomMap,
                agglomProcIDs,
                localMesh,
                otherMeshes,
                g.cellOffsets,
                faceOffsets,
                g.faceMap,
                g.boundaryMap,
                g.boundaryFaceMap
            )
        );
        levelMesh_.set(L_, coarsestMesh_.get());
    }
    else
    {
        levelMesh_.set(L_, nullptr);
    }

    gather_.set(L_, gPtr);
}


void Foam::blockGAMG::allocateWork()
{
    const label nLev = L_ + 1;
    for
    (
        auto* v
      : {&b_, &e_, &r_, &t_, &e2_, &z1_, &q1_, &z2_, &q2_, &rk_, &gbuf_,
         &sd_, &sq_}
    )
    {
        v->clear();
        v->resize(nLev);
    }
    for (label l = 0; l < nLev; ++l)
    {
        if (!hasLevel(l))
        {
            continue;
        }
        const label n = matrixLevel(l).nRows();
        r_[l].resize(n, Zero);
        if (l < L_ && scaleAt(l))
        {
            sd_[l].resize(n, Zero);
            sq_[l].resize(n, Zero);
        }
        if (l > 0)
        {
            b_[l].resize(n, Zero);
            e_[l].resize(n, Zero);
            t_[l].resize(n, Zero);
            e2_[l].resize(n, Zero);
        }
        if (l < L_ && gather_.set(l + 1))
        {
            gbuf_[l].resize(blockDim*aggPtr_->nCells(l), Zero);
        }
    }
    if (cycle_ == cycleKind::K)
    {
        allocateKWork();
    }
}


void Foam::blockGAMG::allocateKWork() const
{
    for (label l = 1; l <= L_; ++l)
    {
        if (hasLevel(l) && z1_[l].empty())
        {
            const label n = matrixLevel(l).nRows();
            z1_[l].resize(n, Zero);
            q1_[l].resize(n, Zero);
            z2_[l].resize(n, Zero);
            q2_[l].resize(n, Zero);
            rk_[l].resize(n, Zero);
        }
    }
}


void Foam::blockGAMG::restrictMatrix(const label fineLevel)
{
    const GAMGAgglomeration& agg = *aggPtr_;
    const blockLduMatrix4& F = matrixLevel(fineLevel);

    // Local coarse sizes (before any processor agglomeration of the
    // coarse level)
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

    // Narrow the local coarse operator
    blockScalarList lDiag(cDiag.size());
    blockScalarList lUpper(cUpper.size());
    blockScalarList lLower(cLower.size());
    forAll(lDiag, i)
    {
        lDiag[i] = narrow(cDiag[i]);
    }
    forAll(lUpper, i)
    {
        lUpper[i] = narrow(cUpper[i]);
        lLower[i] = narrow(cLower[i]);
    }

    // Processor interface coefficients, per fine block interface (the local
    // coarse interface has the same native index)
    const labelListList& patchFaceRestrict =
        agg.patchFaceRestrictAddressing(fineLevel);
    const labelList& nPatchFaces = agg.nPatchFaces(fineLevel);

    const label nFI = F.interfaces().size();
    labelList lIfaceIndex(nFI);
    List<blockScalarList> lIfaceCoeffs(nFI);

    forAll(F.interfaces(), j)
    {
        const label inti = F.interfaces()[j].index();
        lIfaceIndex[j] = inti;

        const blockScalarList& fC = F.interfaceCoeffs(j);
        const labelList& pfr = patchFaceRestrict[inti];

        reduceScalarList sumC(blockSize*nPatchFaces[inti], Zero);
        forAll(pfr, pf)
        {
            const label cpf = pfr[pf];
            for (label k = 0; k < blockSize; ++k)
            {
                sumC[cpf*blockSize + k] += toDouble(fC[pf*blockSize + k]);
            }
        }

        blockScalarList& cC = lIfaceCoeffs[j];
        cC.resize_nocopy(sumC.size());
        forAll(cC, i)
        {
            cC[i] = narrow(sumC[i]);
        }
    }

    if (gather_.set(fineLevel + 1))
    {
        // Processor-agglomerated coarse level (6.3.4)
        gatherMatrix
        (
            fineLevel + 1,
            lDiag,
            lUpper,
            lLower,
            lIfaceIndex,
            lIfaceCoeffs
        );
        return;
    }

    blockLduMatrix4& C = coarse_[fineLevel];

    if (C.nCells() != nCC || C.nFaces() != nCF)
    {
        FatalErrorInFunction
            << "Coarse level " << fineLevel + 1 << " has " << C.nCells()
            << " cells / " << C.nFaces() << " faces, agglomeration "
            << nCC << " / " << nCF << abort(FatalError);
    }

    C.diag() = lDiag;
    C.upper() = lUpper;
    C.lower() = lLower;

    forAll(C.interfaces(), ci)
    {
        const label inti = C.interfaces()[ci].index();
        const label j = lIfaceIndex.find(inti);
        if (j < 0)
        {
            FatalErrorInFunction
                << "No fine-level block interface for coarse interface "
                << inti << " on level " << fineLevel + 1
                << abort(FatalError);
        }
        blockScalarList& cC = C.interfaceCoeffs(ci);
        if (cC.size() != lIfaceCoeffs[j].size())
        {
            FatalErrorInFunction
                << "Coarse interface " << inti << " on level "
                << fineLevel + 1 << " has " << cC.size()/blockSize
                << " faces, agglomeration " << nPatchFaces[inti]
                << abort(FatalError);
        }
        cC = lIfaceCoeffs[j];
    }

    C.markUpdated();
}


void Foam::blockGAMG::gatherMatrix
(
    const label l,
    const blockScalarList& lDiag,
    const blockScalarList& lUpper,
    const blockScalarList& lLower,
    const labelList& lIfaceIndex,
    const List<blockScalarList>& lIfaceCoeffs
)
{
    const gatherLevel& g = gather_[l];

    // Send the local coarse operators to the group master (as native
    // GAMGSolver::gatherMatrices, 16 coefficients per block)
    PstreamBuffers pBufs(g.comm);

    if (!g.master)
    {
        UOPstream toMaster(UPstream::masterNo(), pBufs);
        toMaster
            << lDiag << lUpper << lLower << lIfaceIndex << lIfaceCoeffs;
    }

    pBufs.finishedGathers();

    if (!g.master)
    {
        return;
    }

    // Assemble on the master (as native GAMGSolver::procAgglomerateMatrix)
    blockLduMatrix4& C = coarse_[l - 1];
    C.clear();

    // Native interface index of the gathered level -> block interface
    labelList blockIface(levelMesh_[l].interfaces().size(), -1);
    forAll(C.interfaces(), ci)
    {
        blockIface[C.interfaces()[ci].index()] = ci;
    }

    auto copyBlock = [](blockScalar* dst, const blockScalar* src)
    {
        std::copy(src, src + blockSize, dst);
    };

    const label nProcs = UPstream::nProcs(g.comm);

    for (label proci = 0; proci < nProcs; ++proci)
    {
        blockScalarList rDiag, rUpper, rLower;
        labelList rIfaceIndex;
        List<blockScalarList> rIfaceCoeffs;

        if (proci > 0)
        {
            UIPstream fromProc(proci, pBufs);
            fromProc
                >> rDiag >> rUpper >> rLower >> rIfaceIndex >> rIfaceCoeffs;
        }

        const blockScalarList& pDiag = (proci == 0 ? lDiag : rDiag);
        const blockScalarList& pUpper = (proci == 0 ? lUpper : rUpper);
        const blockScalarList& pLower = (proci == 0 ? lLower : rLower);
        const labelList& pIfaceIndex =
            (proci == 0 ? lIfaceIndex : rIfaceIndex);
        const List<blockScalarList>& pIfaceCoeffs =
            (proci == 0 ? lIfaceCoeffs : rIfaceCoeffs);

        // Diagonal: cells in processor order
        const label nCellsP = g.cellOffsets[proci + 1] - g.cellOffsets[proci];
        const labelList& fMap = g.faceMap[proci];
        if
        (
            pDiag.size() != blockSize*nCellsP
         || pUpper.size() != blockSize*fMap.size()
         || pLower.size() != blockSize*fMap.size()
        )
        {
            FatalErrorInFunction
                << "Level " << l << ": sizes from group rank " << proci
                << " (" << pDiag.size()/blockSize << " cells, "
                << pUpper.size()/blockSize << " faces) do not match the"
                << " agglomeration (" << nCellsP << ", " << fMap.size()
                << ")" << abort(FatalError);
        }
        std::copy
        (
            pDiag.cbegin(),
            pDiag.cend(),
            C.diag().begin() + blockSize*g.cellOffsets[proci]
        );

        // Internal faces (orientation preserved by the gathering; the
        // negative branch is defensive and swaps upper/lower)
        forAll(fMap, facei)
        {
            const label m = fMap[facei];
            const blockScalar* pu = pUpper.cdata() + facei*blockSize;
            const blockScalar* pl = pLower.cdata() + facei*blockSize;
            if (m >= 0)
            {
                copyBlock(C.upper().data() + m*blockSize, pu);
                copyBlock(C.lower().data() + m*blockSize, pl);
            }
            else
            {
                const label a = -m - 1;
                copyBlock(C.upper().data() + a*blockSize, pl);
                copyBlock(C.lower().data() + a*blockSize, pu);
            }
        }

        // Processor interfaces: kept ones map to an interface of the
        // gathered level, merged ones become internal faces. Block
        // coefficients are actual matrix entries A_PN (blockLduInterface):
        // the owner side (map >= 0, its cell is the lower cell) gives the
        // upper block, the neighbour side (map < 0) the lower block.
        const labelList& bMap = g.boundaryMap[proci];
        const labelListList& bfMaps = g.boundaryFaceMap[proci];

        forAll(pIfaceIndex, j)
        {
            const label inti = pIfaceIndex[j];
            if (inti < 0 || inti >= bMap.size() || inti >= bfMaps.size())
            {
                FatalErrorInFunction
                    << "Level " << l << ": interface " << inti
                    << " of group rank " << proci
                    << " not in the boundary map" << abort(FatalError);
            }

            const labelList& bfMap = bfMaps[inti];
            const blockScalarList& pc = pIfaceCoeffs[j];
            if (pc.size() != blockSize*bfMap.size())
            {
                FatalErrorInFunction
                    << "Level " << l << ": interface " << inti
                    << " of group rank " << proci << " has "
                    << pc.size()/blockSize << " faces, boundary face map "
                    << bfMap.size() << abort(FatalError);
            }

            const label allInti = bMap[inti];

            if (allInti >= 0)
            {
                const label ci =
                (
                    allInti < blockIface.size() ? blockIface[allInti] : -1
                );
                if (ci < 0)
                {
                    FatalErrorInFunction
                        << "Level " << l << ": kept interface " << allInti
                        << " (from interface " << inti << " of group rank "
                        << proci << ") is not a block-coupled interface"
                        << abort(FatalError);
                }
                blockScalarList& cc = C.interfaceCoeffs(ci);
                forAll(bfMap, k)
                {
                    copyBlock
                    (
                        cc.data() + bfMap[k]*blockSize,
                        pc.cdata() + k*blockSize
                    );
                }
            }
            else
            {
                forAll(bfMap, k)
                {
                    const label m = bfMap[k];
                    if (m >= 0)
                    {
                        copyBlock
                        (
                            C.upper().data() + m*blockSize,
                            pc.cdata() + k*blockSize
                        );
                    }
                    else
                    {
                        copyBlock
                        (
                            C.lower().data() + (-m - 1)*blockSize,
                            pc.cdata() + k*blockSize
                        );
                    }
                }
            }
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
    const bool gathered = gather_.set(l + 1);

    // Local restriction (into the pre-gather buffer if gathered)
    blockScalarUList& lc = (gathered ? gbuf_[l] : bc);
    lc = Zero;
    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockDim; ++k)
        {
            lc[c*blockDim + k] += r[celli*blockDim + k];
        }
    }

    if (gathered)
    {
        // As native restrictField(..., procAgglom = true), 4 values per cell
        const gatherLevel& g = gather_[l + 1];
        globalIndex::gather
        (
            g.vecOffsets,
            g.comm,
            g.procIDs,
            lc,
            bc,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );
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
    const bool gathered = gather_.set(l + 1);

    if (gathered)
    {
        // As native prolongField(..., procAgglom = true)
        const gatherLevel& g = gather_[l + 1];
        globalIndex::scatter
        (
            g.vecOffsets,
            g.comm,
            g.procIDs,
            e,
            gbuf_[l],
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );
    }

    const blockScalarUList& ec = (gathered ? gbuf_[l] : e);
    forAll(restrictAddr, celli)
    {
        const label c = restrictAddr[celli];
        for (label k = 0; k < blockDim; ++k)
        {
            x[celli*blockDim + k] += ec[c*blockDim + k];
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
        if (diagOn(2))
        {
            diag_->coarseSolve(perf.nIterations, perf.finalResidual, false);
        }
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
    if (diagOn(2))
    {
        diag_->coarseSolve(0, 0, true);
    }
}


void Foam::blockGAMG::smoothLogged
(
    const label l,
    blockScalarUList& x,
    const blockScalarUList& b,
    const label nSweeps,
    DynamicList<doubleScalar>& norms
) const
{
    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        smoothers_[l].smooth(x, b, 1);
        norms.append(residualNorm(l, x, b));
    }
}


Foam::doubleScalar Foam::blockGAMG::residualNorm
(
    const label l,
    const blockScalarUList& x,
    const blockScalarUList& b
) const
{
    const doubleScalar t0 = diagnostics::clock();
    const blockLduMatrix4& A = matrixLevel(l);
    A.residual(r_[l], x, b);
    const doubleScalar nrm = doubleReduce::norm2(r_[l], A.comm());
    diag_->addDiagTime(diagnostics::clock() - t0);
    return nrm;
}


void Foam::blockGAMG::diagOperatorStats() const
{
    // Per scalar row i of the 4x4-block system:
    //   dominance_i = |a_ii| / sum_(j != i) |a_ij|
    // over the diagonal block and the internal-face blocks (processor
    // interface coefficients excluded); rows without off-diagonal entries
    // (e.g. the pressure-reference row) are counted separately.
    for (label l = 0; l <= L_; ++l)
    {
        if (!hasLevel(l) || (l > 0 && !coarse_.set(l - 1)))
        {
            continue;
        }
        const blockLduMatrix4& A = matrixLevel(l);
        const label nC = A.nCells();
        const label nR = blockDim*nC;
        List<doubleScalar> dg(nR, Zero);
        List<doubleScalar> off(nR, Zero);

        const blockScalar* D = A.diag().cdata();
        for (label c = 0; c < nC; ++c)
        {
            for (label r = 0; r < blockDim; ++r)
            {
                for (label k = 0; k < blockDim; ++k)
                {
                    const doubleScalar a =
                        std::abs(toDouble(D[c*blockSize + r*blockDim + k]));
                    if (k == r)
                    {
                        dg[c*blockDim + r] = a;
                    }
                    else
                    {
                        off[c*blockDim + r] += a;
                    }
                }
            }
        }

        const labelUList& lo = A.lduAddr().lowerAddr();
        const labelUList& up = A.lduAddr().upperAddr();
        const blockScalar* Up = A.upper().cdata();
        const blockScalar* Lw = A.lower().cdata();
        forAll(lo, f)
        {
            const label P = lo[f];
            const label N = up[f];
            for (label r = 0; r < blockDim; ++r)
            {
                for (label k = 0; k < blockDim; ++k)
                {
                    const label i = f*blockSize + r*blockDim + k;
                    off[P*blockDim + r] += std::abs(toDouble(Up[i]));
                    off[N*blockDim + r] += std::abs(toDouble(Lw[i]));
                }
            }
        }

        std::vector<doubleScalar> ratio;
        ratio.reserve(std::size_t(nR));
        label nNoOff = 0;
        for (label i = 0; i < nR; ++i)
        {
            if (off[i] > 0)
            {
                ratio.push_back(dg[i]/off[i]);
            }
            else
            {
                ++nNoOff;
            }
        }

        doubleScalar dMin = -1;
        doubleScalar dMed = -1;
        if (!ratio.empty())
        {
            dMin = *std::min_element(ratio.begin(), ratio.end());
            const auto mid = std::ptrdiff_t(ratio.size()/2);
            std::nth_element(ratio.begin(), ratio.begin() + mid, ratio.end());
            dMed = ratio[std::size_t(mid)];
        }
        diag_->operatorStats(l, nC, dMin, dMed, nNoOff);
    }
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

    // Diagnostics (level 2: per-level norms, level 3: per sweep); the
    // lists stay empty (no allocation) below level 3
    const bool d2 = diagOn(2);
    const bool d3 = diagOn(3);
    doubleScalar preBefore = 0;
    doubleScalar preAfter = 0;
    DynamicList<doubleScalar> preSweeps;
    DynamicList<doubleScalar> postSweeps;

    x = Zero;
    if (nPre > 0)
    {
        if (d3)
        {
            smoothLogged(l, x, b, nPre, preSweeps);
        }
        else
        {
            smoothers_[l].smooth(x, b, nPre);
        }
    }

    A.residual(r_[l], x, b);

    if (d2)
    {
        // ||b|| (x = 0 before smoothing) and ||r|| after, one reduction
        const doubleScalar t0 = diagnostics::clock();
        reduceScalar s2[2] =
        {
            doubleReduce::localSumSqr(b),
            doubleReduce::localSumSqr(r_[l])
        };
        doubleReduce::parSum(s2, 2, A.comm());
        preBefore = std::sqrt(s2[0]);
        preAfter = std::sqrt(s2[1]);
        diag_->addDiagTime(diagnostics::clock() - t0);
    }

    // Collective within the group if level l+1 is processor-agglomerated
    restrictVector(l, r_[l], b_[l + 1]);

    blockScalarList& ec = e_[l + 1];

    // Ranks agglomerated away (no level l+1) skip the coarse correction and
    // only take part in the gather above and the scatter below
    if (hasLevel(l + 1))
    {
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
    }

    if (scaleAt(l))
    {
        // Minimal-residual scaling of the prolongated correction d = P e:
        // x += alpha d, alpha = <A d, r>/<A d, A d>, r = b - A x (pre-smoothed)
        blockScalarList& d = sd_[l];
        blockScalarList& q = sq_[l];
        d = Zero;
        prolongAdd(l, ec, d);
        A.Amul(q, d);
        const FixedList<reduceScalar, 2> s =
            doubleReduce::dot2(q, r_[l], q, q, A.comm());
        // GUARD: a vanishing correction is added unscaled (alpha 1)
        const blockScalar alpha =
            narrow(s[1] > doubleScalarVSMALL ? s[0]/s[1] : 1.0);
        const label n = x.size();
        for (label i = 0; i < n; ++i)
        {
            x[i] += alpha*d[i];
        }
    }
    else
    {
        prolongAdd(l, ec, x);
    }

    // r_[l] is free from here on (scratch of the diagnostics norms)
    const doubleScalar postBefore = (d2 ? residualNorm(l, x, b) : 0);

    if (nPost > 0)
    {
        if (d3)
        {
            smoothLogged(l, x, b, nPost, postSweeps);
        }
        else
        {
            smoothers_[l].smooth(x, b, nPost);
        }
    }

    if (d2)
    {
        const doubleScalar postAfter =
        (
            postSweeps.size()
          ? postSweeps.last()
          : (nPost > 0 ? residualNorm(l, x, b) : postBefore)
        );
        const doubleScalar t0 = diagnostics::clock();
        diag_->levelVisit
        (
            l, preBefore, preAfter, postBefore, postAfter,
            preSweeps, postSweeps
        );
        diag_->addDiagTime(diagnostics::clock() - t0);
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

    // Only ranks holding level l get here; all reductions use its
    // communicator
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

    // Diagnostics: r1norm -1 = not computed (single-step K cycle)
    if (kMaxSteps_ < 2)
    {
        if (diagOn(2))
        {
            diag_->kStep(l, r0norm, -1, kThreshold_*r0norm, false, a1, 0);
        }
        return;
    }
    const reduceScalar r1norm = doubleReduce::norm2(r1, comm);
    if (r1norm <= kThreshold_*r0norm)
    {
        if (diagOn(2))
        {
            diag_->kStep
            (
                l, r0norm, r1norm, kThreshold_*r0norm, false, a1, 0
            );
        }
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

    if (diagOn(2))
    {
        diag_->kStep(l, r0norm, r1norm, kThreshold_*r0norm, true, a1, a2);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::blockGAMG::nSingularDiag() const
{
    label n = 0;
    forAll(smoothers_, l)
    {
        if (smoothers_.set(l))
        {
            n += smoothers_[l].nSingularDiag();
        }
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
    if (c == cycleKind::K)
    {
        // K workspace may not have been allocated for the start-up cycle
        allocateKWork();
    }
    cycle_ = c;
    return true;
}


void Foam::blockGAMG::update()
{
    ++nUpdates_;
    const bool matrixWeights = (aggWeights_ != "geometric");
    if
    (
        matrixWeights
     && (
            !aggFromMatrix_
         || (reaggInterval_ > 0 && nUpdates_ % reaggInterval_ == 0)
        )
    )
    {
        // First update (the matrix exists only now) or re-agglomeration:
        // pair agglomeration on the block-matrix weights
        clearHierarchy();
        aggPtr_ = nullptr;
        aggFromMatrix_ = true;
        buildHierarchy();
    }
    else if (!dict_.get<bool>("cacheAgglomeration"))
    {
        // Rebuild everything on a fresh agglomeration (the matrices and
        // smoothers reference the agglomeration's addressing)
        const label nLevOld = nLevels();
        clearHierarchy();
        fine_.mesh().thisDb().checkOut
        (
            const_cast<GAMGAgglomeration*>(aggPtr_)
        );
        aggPtr_ = nullptr;
        buildHierarchy();
        if (nLevels() != nLevOld)
        {
            FatalErrorInFunction
                << "Level count changed on re-agglomeration of a static mesh"
                << abort(FatalError);
        }
    }

    // Levels in order: a group master assembles level l+1 from the group
    // before restricting it further
    for (label l = 0; l < L_; ++l)
    {
        if (hasLevel(l))
        {
            restrictMatrix(l);
        }
    }

    forAll(smoothers_, l)
    {
        if (smoothers_.set(l))
        {
            smoothers_[l].update();
        }
    }

    factoriseDense();

    if (diagOn(3))
    {
        const doubleScalar t0 = diagnostics::clock();
        diagOperatorStats();
        diag_->addDiagTime(diagnostics::clock() - t0);
    }
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
    os  << "blockGAMG: levels " << nLevels()
        << ", mergeLevels " << mergeLevelsUsed_
        << ", C_op " << Cop_
        << ", cycle " << cycleName(cycle_)
        << ", processorAgglomerator " << procAgglomType_ << nl
        << "blockGAMG: cells per level " << flatOutput(cellsPerLevel_) << nl
        << "blockGAMG: ranks per level " << flatOutput(ranksPerLevel_) << nl
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
    d.add("cellsPerLevel", cellsPerLevel_);
    d.add("ranksPerLevel", ranksPerLevel_);
    d.add("processorAgglomerator", procAgglomType_);
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
    dict.set("processorAgglomerator", procAgglomType_);
}


// ************************************************************************* //
