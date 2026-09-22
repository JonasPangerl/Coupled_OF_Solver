/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::diagnosticFields

Description
    Write-time diagnostic fields (amendment C6, spec 12 addition). At every
    write time and on the 9.3 abort path (into <n>_lastValid):
      localDt       dt_P after all caps (PTC CFL, remediation and zonal
                    factors, local limit 7.3) [s]
      localCFL      dt_P (1/2 sum_f |phi_f|)/V_P with the phi of the written
                    state [-]
      cflFactorEff  effective per-cell CFL multiplier (static, dynamic and
                    zonal factors combined; product of the zonal factors) [-]
      betaEff       effective per-cell convection blending factor (global
                    start-up value times the static/dynamic/zonal rules) [-]
      USFD          the SFD filtered velocity, only if a field was handed in
                    with setSFD() (C3 hook)
    remediationFlag (bits 1/2/4/8/16) is written by remediation::write at
    the same places.

    dt_P, cflFactorEff and betaEff are those of the last assembled trial
    step (record()); on the abort path that is the aborted step.

    The fields use extrapolatedCalculated boundaries (patch value = cell
    value) and are written through the file handler, so that a non-time
    instance such as <n>_lastValid is honoured (regIOobject::writeObject
    redirects it to the current time name in v2606).

    Dictionary: coupled.diagnosticFields { enabled yes; }

SourceFiles
    diagnosticFields.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_diagnosticFields_H
#define coupledFoam_diagnosticFields_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class diagnosticFields
{
    // Private Data

        const fvMesh& mesh_;

        bool enabled_;

        //- V_P/dt_P of the last assembled trial
        scalarField rDeltaTV_;

        //- Effective CFL and beta factors of the last assembled trial
        scalarField cflFactor_;
        scalarField beta_;

        //- SFD filtered velocity (C3 hook; nullptr: not written)
        const volVectorField* USFD_;


    // Private Member Functions

        //- A dimensionless or dimensioned cell field for the instance
        tmp<volScalarField> makeField
        (
            const word& name,
            const word& instance,
            const dimensionSet& dims,
            const scalarField& values
        ) const;

        //- Write through the file handler into io.instance()
        static void writeInstance(const regIOobject& io);


public:

    // Constructors

        diagnosticFields(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        bool enabled() const noexcept { return enabled_; }

        //- Store the per-cell data of the current trial step (cheap copies)
        void record
        (
            const scalarField& rDeltaTV,
            const scalarField& cflFactor,
            const scalarField& beta
        );

        //- C3 hook: the SFD filtered velocity to be written as USFD while
        //  SFD is active (nullptr: none)
        void setSFD(const volVectorField* Ubar) noexcept { USFD_ = Ubar; }

        //- Write the fields into the given instance (collective)
        void write(const surfaceScalarField& phi, const word& instance) const;

        //- Write the fields into the current time
        void write(const surfaceScalarField& phi) const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
