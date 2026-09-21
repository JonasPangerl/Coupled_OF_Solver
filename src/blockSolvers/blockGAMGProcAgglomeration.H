/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGAMGProcAgglomeration

Description
    Native GAMGProcAgglomeration (runtime name "blockGAMGRule") that
    implements the processor-agglomeration rule of spec 6.3.4:

        T(N) = nRanks                        if N >= cellsPerRank*nRanks
             = max(1, nRanks/divisor)        if cellsPerRank <= N
                                                < cellsPerRank*nRanks
             = 1                             if N < cellsPerRank

    where N is the global cell count of a level (measured on the finished
    local agglomeration) and nRanks the rank count of the finest mesh.
    Defaults cellsPerRank = coupledDefaults::procAgglomCellsPerRank (5000),
    divisor = coupledDefaults::procAgglomDivisor (4); dictionary keys
    procAgglomCellsPerRank, procAgglomDivisor.

    It works exactly like the native "manual" agglomerator (same call of
    GAMGProcAgglomeration::agglomerate, i.e. native mesh gathering by
    procAgglomerateLduAddressing and re-agglomeration of the coarser levels
    on the masters), but the per-level processor maps are computed from the
    rule instead of being read from the dictionary: at every level l in
    1..L-1 whose target T(N_l) is below the current rank count R, the R
    ranks of the level communicator are grouped into T contiguous, balanced
    blocks (rank i -> group (i*T)/R, master = lowest rank of a group).

    Limits of the native mechanism (handled by blockGAMG):
      - level 0 (the decomposed mesh) is never agglomerated;
      - the coarsest level L cannot be gathered natively (the native
        gathering also gathers the restriction addressing below the
        level); if the rule requires fewer ranks on L than on L-1,
        blockGAMG gathers L itself with the same native lduPrimitiveMesh
        gathering (see blockGAMG::gatherCoarsest).

    The class is implemented in blockGAMG.C (not in a file of its own), so
    that src/Make/files does not have to change. Moving it into
    blockGAMGProcAgglomeration.C is a pure cut-and-paste plus one line in
    src/Make/files.

SourceFiles
    blockGAMG.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGAMGProcAgglomeration_H
#define coupledFoam_blockGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"
#include "blockScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        //- Cells per rank threshold (and single-rank threshold)
        const label cellsPerRank_;

        //- Rank divisor of the intermediate band
        const label divisor_;


public:

    //- Runtime type information
    TypeName("blockGAMGRule");


    // Constructors

        //- Construct given agglomerator and controls
        blockGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );

        //- No copy
        blockGAMGProcAgglomeration(const blockGAMGProcAgglomeration&) = delete;
        void operator=(const blockGAMGProcAgglomeration&) = delete;


    //- Destructor
    virtual ~blockGAMGProcAgglomeration() = default;


    // Static Functions

        //- Target rank count of a level with nCells global cells (rule 6.3.4)
        static label targetRanks
        (
            const reduceScalar nCells,
            const label nRanks,
            const label cellsPerRank,
            const label divisor
        );

        //- Processor map of nProcs ranks onto nGroups contiguous balanced
        //  groups: rank i -> (i*nGroups)/nProcs
        static labelList groupMap(const label nProcs, const label nGroups);


    // Member Functions

        //- Apply the rule to levels 1..L-1. Returns true.
        virtual bool agglomerate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
