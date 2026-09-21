/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::remediation

Description
    Two-tier remediation cell sets (spec 8).

    Static set (8.1), built once from native per-face mesh-quality fields
    (polyMeshTools::faceOrthogonality/faceSkewness/volRatio,
    primitiveMeshTools::cellClosedness for the aspect ratio). A cell is
    marked if any of its faces exceeds nonOrthThreshold (deg) or
    skewThreshold, the volume ratio to any neighbour exceeds
    volRatioThreshold, or its aspect ratio exceeds aspectThreshold.
    Treatment: beta = static.beta, non-orthogonal limiter static value
    (via coupledAssembler::setStaticCells), CFL factor static.cflFactor.
    The "gradient limiter forced to cellLimited 1" of 8.1 is recorded in
    DECISIONS.md D-018.

    Dynamic set (8.2), updated every iteration. Marked if any of
        |U_P| > cU Uref,   |p_P| > cp pref,
        |U_P - Ubar_N(P)| > cSpike Uref  (face-neighbour mean),
        non-finite U, p, k, omega,
    or if passed in by the line search / sentinel. Marked cells plus nLayers
    layers of face neighbours (across processor faces as well) get beta = 0,
    CFL factor dynamic.cflFactor and optional increment clipping
    |U_P + dU_P - Ubar_N| <= cSpike Uref (the field itself is never
    overwritten). Hysteresis: a cell leaves the set after nQuietIters
    consecutive iterations without being re-marked.

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

    Output: cellSets remediationStatic, remediationDynamic and the
    volScalarField remediationFlag (0 none, 1 static, 2 dynamic, 3 both);
    with zonal factors enabled additionally the volScalarField
    zonalCflFactor.

SourceFiles
    remediation.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_remediation_H
#define coupledFoam_remediation_H

#include "fvMesh.H"
#include "volFields.H"
#include "blockScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class remediation
{
    // Private Data

        const fvMesh& mesh_;

        // Static controls
        bool staticEnabled_;
        scalar nonOrthThreshold_;
        scalar skewThreshold_;
        scalar volRatioThreshold_;
        scalar aspectThreshold_;
        scalar staticBeta_;
        scalar staticCflFactor_;

        // Dynamic controls
        bool dynamicEnabled_;
        scalar cU_;
        scalar cp_;
        scalar cSpike_;
        label nLayers_;
        label nQuietIters_;
        scalar dynamicCflFactor_;
        bool clipToNeighbourMean_;
        scalar warnFraction_;
        label warnInterval_;

        // State

            //- Static membership
            boolList isStatic_;

            //- Dynamic age per cell: -1 not in the set, else iterations
            //  since the cell was last marked
            labelList age_;

            //- Cells marked externally since the last update
            //  (line search, sentinel)
            boolList pending_;

            //- Global counts
            label nStatic_;
            label nDynamic_;

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


public:

    // Constructors

        remediation(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        //- Build the static set (8.1). Prints a summary.
        void buildStatic();

        //- Mark cells for the dynamic set at the next update
        void markDynamic(const labelUList& cells);

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
        //  (start-up switch): static cells static.beta, dynamic cells 0,
        //  multiplied by the zonal beta factors
        tmp<scalarField> beta(const scalar betaGlobal) const;

        //- Per-cell CFL factor (1, static, dynamic; the smallest applies),
        //  multiplied by the zonal CFL factors
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
        label nStatic() const noexcept { return nStatic_; }
        label nDynamic() const noexcept { return nDynamic_; }
        bool zonalEnabled() const noexcept { return zonalEnabled_; }

        //- Write remediationFlag and the two cellSets for the current time
        //  (and zonalCflFactor if zonal factors are enabled)
        void write() const;

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
