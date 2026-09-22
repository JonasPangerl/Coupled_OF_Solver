/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockSmootherPrecon

Description
    Single-level preconditioner (precond-research): nSweeps sweeps of a
    block smoother (keyword smoother, as in blockGAMG) from a zero guess,
        preconditioner blockSmoother; smoother blockILU0; nSweeps 1;
    With blockILU0 and nSweeps 1 this is the block ILU0 preconditioner.
    Used to separate the smoother from the coarse-grid correction.

SourceFiles
    blockSmootherPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockSmootherPrecon_H
#define coupledFoam_blockSmootherPrecon_H

#include "blockPreconditioner.H"
#include "blockSmoother.H"

namespace Foam
{

class blockSmootherPrecon
:
    public blockPreconditioner
{
    autoPtr<blockSmoother> smoother_;
    label nSweeps_;

public:

    TypeName("blockSmoother");

    blockSmootherPrecon(const blockSolver& solver, const dictionary& dict);

    virtual ~blockSmootherPrecon() = default;

    virtual void update();

    virtual void precondition
    (
        blockScalarUList& w,
        const blockScalarUList& r
    ) const;

    virtual label nSingularDiag() const
    {
        return smoother_->nSingularDiag();
    }

    virtual void writeSettings(dictionary& dict) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
