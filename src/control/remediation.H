/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::remediation

Description
    Two-tier remediation cell sets (spec 8).

    Static set (8.1): cells selected in advance, built once at start-up and
    split into four categories with their own switches, criteria and
    treatment (D-066; the pre-D-066 remediation.static entries map onto
    them, see below):
    \verbatim
    meshQuality     // 8.1 quality criteria (full treatment)
    {
        enabled yes; nonOrthThreshold 85; skewThreshold 6;
        volRatioThreshold 30; aspectThreshold 2000;
        beta 0; cflFactor 0.5; gradLimiter yes; nonOrthLimiter yes;
    }
    badMesh         // severe defects (full treatment)
    {
        enabled yes; volumeJump yes; volJumpThreshold 0.98;
        closednessThreshold 1e-6;       // 0 = off
        nonOrthThreshold 0; skewThreshold 0;    // severe values, 0 = off
        beta 0; cflFactor 0.5; gradLimiter yes; nonOrthLimiter yes;
    }
    processor       // cells at processor patches (mild treatment)
    {
        enabled yes; procAMI yes; nLayers 0;
        beta 1; cflFactor 0.5; gradLimiter no; nonOrthLimiter no;
    }
    wall            // cells at wall patches (mild treatment)
    {
        enabled yes; wallStarved no; maxWallInternalFaces 2;
        patches (); nLayers 1;          // wordRes, layer 1 = patch cells
        beta 1; cflFactor 0.5; gradLimiter no; nonOrthLimiter no;
    }
    limitedNonOrthCoeff 0.2;            // nonOrthLimiter cells
    limitedGradScheme "cellLimited Gauss linear 1";     // gradLimiter cells
    \endverbatim
    meshQuality uses the native per-face mesh-quality fields
    (polyMeshTools::faceOrthogonality/faceSkewness/volRatio,
    primitiveMeshTools::cellClosedness for the aspect ratio): a cell is
    marked if any of its faces exceeds nonOrthThreshold (deg) or
    skewThreshold, the volume ratio to any neighbour exceeds
    volRatioThreshold, or its aspect ratio exceeds aspectThreshold.
    badMesh: C1 volumeJump, cell openness (cellClosedness) above
    closednessThreshold, faces above the severe nonOrth/skew thresholds.
    processor: C1 procAMI and nLayers face layers from the processor
    patches (decomposition dependent). wall: C1 wallStarved (default off,
    D-061) and nLayers face layers from the patches of wall.patches
    (wordRes: names, regular expressions, groups such as wall).

    Treatment: beta (upper bound of the convection blending factor),
    cflFactor (local time step factor), gradLimiter (the p gradient of the
    Rhie-Chow terms from limitedGradScheme, D-018), nonOrthLimiter (the
    faces of the cell use limitedNonOrthCoeff instead of the global
    nonOrthLimiter; via coupledAssembler::setLimitedCells). A cell in
    several categories gets the most restrictive value of each: min beta,
    min cflFactor, OR of the limiter switches.

    Legacy keywords (pre-D-066, one deprecation note in the log):
    remediation.static { enabled -> all four categories; nonOrthThreshold,
    skewThreshold, volRatioThreshold, aspectThreshold -> meshQuality;
    beta, cflFactor -> meshQuality and badMesh; nonOrthLimiter ->
    limitedNonOrthCoeff; volumeJump, volJumpThreshold -> badMesh; procAMI ->
    processor; wallStarved, maxWallInternalFaces -> wall (mild treatment,
    D-061) }. A keyword in the category dictionary takes precedence.

    Dynamic set (8.2), updated every iteration. Marked if any of
        |U_P| > cU Uref,   |p_P| > cp pref,
        |U_P - Ubar_N(P)| > cSpike Uref  (face-neighbour mean),
        non-finite U, p, k, omega,
    or if passed in by the line search / sentinel. Marked cells plus nLayers
    layers of face neighbours (across processor faces as well) get beta = 0,
    CFL factor dynamic.cflFactor and optional increment clipping
    |U_P + dU_P - Ubar_N| <= cSpike Uref (the field itself is never
    overwritten). Hysteresis: a cell leaves the set after nQuietIters
    consecutive iterations without being re-marked. A cell that entered the
    set stickyAfter times (0 = never, D-055) stays in it for the rest of
    the run; the per-cell entry counters are part of the restart state.
    Release ramp (dynamic.releaseIters N > 0, D-055): a released cell is
    out of the set (no clipping, not counted in nDyn) but for N further
    iterations, with j = N, N-1, ..., 1 iterations of the ramp left, its
    beta is scaled by 1 - j/(N+1) and its CFL factor is
    cflFactor^(j/(N+1)). A cell
    re-marked during the ramp re-enters the set. 0 = the instant release
    of 8.2 (bit-identical).

    Zonal factors (8.3, amendment B6), built once at construction:
    \verbatim
    zonal
    {
        enabled no;
        zones ( { cellZone wheelHousing; cflFactor 0.5; beta 0.5; } );
        patchDistance
        ( { patches (body); nLayers 3; cflFactor 0.7; beta 1.0; } );
    }
    \endverbatim
    Multipliers on the local time step (cflFactor) and on the convection
    blending factor (beta) for cells in the named cellZones or within
    nLayers face-neighbour layers of the named patches (layer 1 = cells
    adjacent to the patch faces, layer k+1 = face neighbours of layer k,
    across processor faces as well). Zone and patch names are wordRe
    (regular expressions and groups allowed). Several matching entries
    multiply. The multipliers apply on top of the static/dynamic rules
    (dynamic cells keep beta = 0). Unknown names are a FatalIOError.

    The dynamic set always overrides the static treatment (beta 0,
    dynamic.cflFactor); the zonal factors multiply on top of both.

    Amendment C1: topological static criteria wallStarved, procAMI and
    volumeJump (staticCriteria.H), sorted into the categories; C5: built-in
    named sets _wallCells ... _remediationStatic in zonal.zones
    (zonalSets.H), plus one per category (_remediationMeshQuality, ...).

    Output: cellSets remediationStatic (union of the categories),
    remediationDynamic, remediationMeshQuality, remediationBadMesh,
    remediationProcessor, remediationWall and the volScalarField
    remediationFlag (bits OR-ed, staticCriteria::flagBits: 1 static,
    2 dynamic, 4 wallStarved, 8 procAMI, 16 volumeJump, 32 meshQuality,
    64 badMesh, 128 processor, 256 wall, 512 closedness, 1024 severe
    quality, 2048 processor layer, 4096 wall layer); with zonal factors
    enabled additionally the volScalarField zonalCflFactor.

SourceFiles
    remediation.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_remediation_H
#define coupledFoam_remediation_H

#include "fvMesh.H"
#include "volFields.H"
#include "blockScalar.H"
#include "FixedList.H"
#include "wordRes.H"
#include "staticCriteria.H"
#include "zonalSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class remediation
{
public:

    //- Pre-selected cell categories (D-066)
    enum category : label
    {
        meshQuality = 0,
        badMesh,
        processor,
        wall,
        nCategories
    };

    //- Switch and treatment of one category
    struct categorySettings
    {
        bool enabled;
        scalar beta;
        scalar cflFactor;
        bool gradLimiter;
        bool nonOrthLimiter;
    };

    //- Category names (dictionary keywords) and remediationFlag bits
    static const char* categoryName(const label cat);
    static label categoryBit(const label cat);


private:

    // Private Data

        const fvMesh& mesh_;

        // Static categories (D-066)
        FixedList<categorySettings, nCategories> cat_;

        // meshQuality criteria (8.1)
        scalar nonOrthThreshold_;
        scalar skewThreshold_;
        scalar volRatioThreshold_;
        scalar aspectThreshold_;

        // badMesh criteria besides volumeJump (0 = off)
        scalar closednessThreshold_;
        scalar severeNonOrthThreshold_;
        scalar severeSkewThreshold_;

        // processor: face layers from the processor patches (0 = none)
        label procNLayers_;

        // wall: face layers from the wall.patches patches
        wordRes wallPatches_;
        label wallNLayers_;

        //- Topological criteria (amendment C1), switches of disabled
        //  categories cleared
        staticCriteria::settings criteria_;

        //- Per-cell static bits (staticCriteria::flagBits without
        //  bitDynamic: criteria, categories and bitStatic)
        labelList staticBits_;

        //- Per-cell treatment of the static cells (meaningful where
        //  isStatic_): min beta, min cflFactor over the cell's categories
        scalarField staticBeta_;
        scalarField staticCfl_;

        //- Per-cell limiter switches (OR over the cell's categories)
        boolList gradLimited_;
        boolList nonOrthLimited_;

        //- Global counts per category
        FixedList<label, nCategories> nCat_;

        // Dynamic controls
        bool dynamicEnabled_;
        scalar cU_;
        scalar cp_;
        scalar cSpike_;
        label nLayers_;
        label nQuietIters_;
        label stickyAfter_;
        label releaseIters_;
        scalar dynamicCflFactor_;
        bool clipToNeighbourMean_;
        scalar warnFraction_;
        label warnInterval_;

        // State

            //- Static membership (union of the categories)
            boolList isStatic_;

            //- Dynamic age per cell: -1 not in the set, else iterations
            //  since the cell was last marked
            labelList age_;

            //- Dynamic entries per cell (transitions into the set); a cell
            //  with stickyAfter_ or more entries is never released
            labelList entries_;

            //- Global number of sticky cells in the set
            label nSticky_;

            //- Release ramp iterations left per cell (0: none)
            labelList rampLeft_;

            //- Global number of cells in a release ramp
            label nRamping_;

            //- Cells marked externally since the last update
            //  (line search, sentinel)
            boolList pending_;

            //- Global counts
            label nStatic_;
            label nDynamic_;

            //- Incremented whenever the dynamic-set membership changes
            //  (global; the Anderson history is flushed on a change)
            label dynamicVersion_;

        // Zonal factors (8.3)

            bool zonalEnabled_;

            //- Per-cell product of the zonal CFL multipliers
            //  (empty if zonal factors are disabled)
            scalarField zonalCfl_;

            //- Per-cell product of the zonal beta multipliers
            //  (empty if zonal factors are disabled)
            scalarField zonalBeta_;

            //- Effective zonal settings (defaults filled in)
            dictionary zonalSettings_;

            //- Built-in entries applied after the static set is built
            //  (_remediationStatic, amendment C5)
            DynamicList<zonalSets::deferredEntry> zonalDeferred_;


    // Private Member Functions

        //- Read the zonal dictionary, build the per-cell multipliers and
        //  print a one-line summary per entry (collective)
        void buildZonal(const dictionary& zonalDict);

        //- Multiply one entry's factors into the per-cell multipliers.
        //  Returns the global number of affected cells.
        label applyZonal
        (
            const boolList& mark,
            const scalar cflFactor,
            const scalar beta
        );

        //- Apply the deferred built-in zonal entries (C5); collective
        void applyDeferredZonal();

        //- Read and check the cflFactor and beta of one zonal entry
        void readZonalFactors
        (
            const dictionary& entryDict,
            const dictionary& zonalDict,
            scalar& cflFactor,
            scalar& beta
        ) const;

        //- Face-neighbour mean of U per cell (processor faces included)
        tmp<vectorField> neighbourMean(const volVectorField& U) const;

        //- Grow a marker by one face layer (processor faces included)
        void growLayer(boolList& mark) const;

        //- Cells within nLayers face layers of the (non-processor) patches
        //  matching patchNames (layer 1 = cells with a face on them).
        //  Unknown names are a FatalIOError naming the context. Collective.
        boolList patchLayers
        (
            const wordRes& patchNames,
            const label nLayers,
            const dictionary& errDict,
            const string& context
        ) const;

        //- Read the categories (with the legacy remediation.static
        //  mapping) from the remediation dictionary
        void readCategories(const dictionary& remDict);


public:

    // Constructors

        remediation(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        //- Build the static set (8.1, categories D-066). Prints a summary.
        //  Collective.
        void buildStatic();

        //- Mark cells for the dynamic set at the next update
        void markDynamic(const labelUList& cells);

        //- Put the marked cells into the dynamic set NOW, with nLayers +
        //  extraLayers halo layers, without ageing any other cell (no step
        //  was accepted). For the sentinel rollback (D-076): the retry has
        //  to be assembled with the offending cells already treated, and
        //  updateDynamic only runs after an accepted step, i.e. never in a
        //  rollback cascade. Independent of startupReference exclude, which
        //  concerns the field classification, not cells the sentinel has
        //  caught. Returns the number of cells that entered (global).
        //  Collective.
        label activatePending(const label extraLayers = 0);

        //- Dynamic update (8.2) at the end of an outer iteration
        void updateDynamic
        (
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField* kPtr,
            const volScalarField* omegaPtr,
            const scalar Uref,
            const scalar pref,
            const label iter
        );

        //- Per-cell convection blending factor given the global value
        //  (start-up switch): static cells at most their category beta,
        //  dynamic cells 0, multiplied by the zonal beta factors
        tmp<scalarField> beta(const scalar betaGlobal) const;

        //- Per-cell CFL factor (1, category, dynamic; the smallest
        //  applies), multiplied by the zonal CFL factors
        tmp<scalarField> cflFactor() const;

        //- Clip the velocity increment in dynamic cells (8.2)
        //  Returns the number of clipped cells (this rank).
        label clipIncrement
        (
            blockScalarUList& dx,
            const volVectorField& U,
            const scalar omega,
            const scalar Uref
        ) const;

        const boolList& isStatic() const noexcept { return isStatic_; }

        //- Cells whose p gradient in the Rhie-Chow terms is limited
        const boolList& gradLimited() const noexcept { return gradLimited_; }

        //- Cells whose faces use the limited non-orthogonal coefficient
        const boolList& nonOrthLimited() const noexcept
        {
            return nonOrthLimited_;
        }

        //- Global number of static cells (union of the categories)
        label nStatic() const noexcept { return nStatic_; }

        //- Global number of cells in a category
        label nCategory(const label cat) const { return nCat_[cat]; }

        label nDynamic() const noexcept { return nDynamic_; }
        label nSticky() const noexcept { return nSticky_; }
        label nRamping() const noexcept { return nRamping_; }

        //- Dynamic age per cell (-1: not in the set) and entry counts
        const labelList& age() const noexcept { return age_; }
        const labelList& entries() const noexcept { return entries_; }
        label dynamicVersion() const noexcept { return dynamicVersion_; }
        bool zonalEnabled() const noexcept { return zonalEnabled_; }

        //- Write remediationFlag and the cellSets (static, dynamic and one
        //  per category) for the current time (and zonalCflFactor if zonal
        //  factors are enabled)
        void write() const;

        //- The same into the given instance (e.g. <iter>_lastValid of the
        //  9.3 abort path, so that no incomplete time directory is created)
        void write(const word& instance) const;

        //- Restart state: dynamic set and ages (local labels)
        void writeState(dictionary& dict) const;
        void readState(const dictionary& dict);

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
