/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::MRFCoupling

Description
    Implicit MRF (spec 5.6, DECISIONS.md D-006).

    Native MRF (absolute velocity, relative flux) adds the frame acceleration
    Omega x U_abs explicitly (MRFZoneList::DDt). Here the same term is made
    implicit: V_P [Omega]_x is added to the 3x3 momentum part of the
    diagonal block, where [Omega]_x u = Omega x u.

    The per-cell matrix is extracted from the native MRFZoneList::DDt applied
    to the three unit vector fields, so any number of zones and a
    time-dependent Omega are handled through the public API only. It is
    re-extracted every outer iteration (cheap: three explicit evaluations).

    Flux convention as native simpleFoam: phi is stored relative
    (MRF.makeRelative after the flux update).

SourceFiles
    MRFCoupling.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_MRFCoupling_H
#define coupledFoam_MRFCoupling_H

#include "fvMesh.H"
#include "IOMRFZoneList.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class MRFCoupling
{
    // Private Data

        const fvMesh& mesh_;

        const IOMRFZoneList& MRF_;

        //- V_P [Omega]_x per cell (zero outside the zones)
        tensorField coeffs_;

        //- True if any MRF zone is active
        bool active_;


public:

    // Constructors

        MRFCoupling(const fvMesh& mesh, const IOMRFZoneList& MRF);


    // Member Functions

        //- Active?
        bool active() const noexcept
        {
            return active_;
        }

        //- The MRF zone list
        const IOMRFZoneList& MRF() const noexcept
        {
            return MRF_;
        }

        //- Re-extract V [Omega]_x (call once per outer iteration)
        void update();

        //- V_P [Omega]_x (empty if not active)
        const tensorField& coeffs() const noexcept
        {
            return coeffs_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
