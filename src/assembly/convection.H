/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::convection

Description
    Momentum convection with deferred correction (spec 5.3a).

    Implicit part: native fvm::div(phi, U) with the low-order scheme of the
    fvSchemes entry "div(phi,U)_upwind" (default in the case templates:
    bounded Gauss upwind), i.e. owner diag += max(phi,0), upper = -max(-phi,0),
    lower = -max(phi,0), neighbour diag += max(-phi,0).

    Explicit part, per unit volume, to be subtracted from the right-hand
    side after multiplication with V_P:
        beta_P * [ fvc::div(phi,U,"div(phi,U)") - fvc::div(phi,U,"div(phi,U)_upwind") ]
    with the per-cell blending factor beta in [0,1] (start-up switch 7.4,
    remediation 8).

SourceFiles
    convection.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_convection_H
#define coupledFoam_convection_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class convection
{
    // Private Data

        const fvMesh& mesh_;

        //- fvSchemes key of the high-order scheme
        const word hoName_;

        //- fvSchemes key of the low-order (implicit) scheme
        const word udName_;


public:

    // Constructors

        explicit convection(const fvMesh& mesh);


    // Member Functions

        const word& hoName() const noexcept { return hoName_; }
        const word& udName() const noexcept { return udName_; }

        //- Implicit low-order operator
        tmp<fvVectorMatrix> implicitPart
        (
            const surfaceScalarField& phi,
            const volVectorField& U
        ) const;

        //- Deferred correction per unit volume, beta-weighted:
        //  beta*(div_HO - div_UD). Zero where beta == 0; the HO evaluation is
        //  skipped entirely if beta is zero everywhere (start-up phase).
        tmp<vectorField> deferredCorrection
        (
            const surfaceScalarField& phi,
            const volVectorField& U,
            const scalarField& beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
