/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockSimplePrecon

Description
    SIMPLE-type block-factorisation preconditioner (precond-research,
    selection name "blockSimple") for the saddle-point structure
        [A G] [u]   [r_u]
        [B C] [p] = [r_p]
    with the approximate factorisation
        [A 0] [I  L^-1 G]
        [B S] [0  I     ],   L = diag(A) (u-u block diagonal entries),
    and the Schur approximation S = schurScale*C (C: the p-p block, the
    compact Rhie-Chow Laplacian with D = V/abar; B L^-1 G is a Laplacian
    of about the same weight, hence the default schurScale 2).

    A and S are solved approximately by one blockGAMG cycle on the
    DECOUPLED block matrix Md = [A 0; 0 S] (u-u and p-p couplings only),
    which is block diagonal in (u, p), so smoothing and coarse operators
    act on a momentum and a pressure-Laplacian problem.

    simpleMode
      sequential  y_u = GAMG(r_u, 0); y_p = GAMG(0, r_p - B y_u)  (2 cycles)
      single      t = L^-1 r_u; (y_u, y_p) = GAMG(r_u, r_p - B t) (1 cycle)
    then u = y_u - L^-1 G y_p, p = y_p.

    Keywords (solver dictionary): simpleMode (single), schurScale (2),
    blockGAMG { ... } (hierarchy on Md). Variable preconditioner: requires
    blockFGMRES.

SourceFiles
    blockSimplePrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockSimplePrecon_H
#define coupledFoam_blockSimplePrecon_H

#include "blockPreconditioner.H"
#include "blockGAMG.H"

namespace Foam
{

class blockSimplePrecon
:
    public blockPreconditioner
{
    //- Decoupled matrix [A 0; 0 S]
    blockLduMatrix4 Md_;

    //- Hierarchy on Md
    blockGAMG gamg_;

    //- Schur scale and mode
    reduceScalar schurScale_;
    bool sequential_;

    //- 1/diag of the u-u block per cell and u component (p entry 0)
    blockScalarList rL_;

    //- Work vectors
    mutable blockScalarList t_, v_, y_, y2_;

public:

    TypeName("blockSimple");

    blockSimplePrecon(const blockSolver& solver, const dictionary& dict);

    virtual ~blockSimplePrecon() = default;

    const blockGAMG& gamg() const noexcept { return gamg_; }

    virtual void update();

    virtual void precondition
    (
        blockScalarUList& w,
        const blockScalarUList& r
    ) const;

    virtual label nSingularDiag() const
    {
        return gamg_.nSingularDiag();
    }

    virtual void writeSettings(dictionary& dict) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
