/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockSmootherPrecon.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(blockSmootherPrecon, 0);
    addToRunTimeSelectionTable
    (
        blockPreconditioner,
        blockSmootherPrecon,
        dictionary
    );
}


Foam::blockSmootherPrecon::blockSmootherPrecon
(
    const blockSolver& solver,
    const dictionary& dict
)
:
    blockPreconditioner(solver),
    smoother_(),
    nSweeps_
    (
        dict.getOrDefault<label>
        (
            "nSweeps",
            coupledDefaults::smootherPreconSweeps
        )
    )
{
    dictionary sd(dict);
    if (!sd.found("smoother"))
    {
        sd.add("smoother", word(coupledDefaults::smootherPreconSmoother));
    }
    smoother_ = blockSmoother::New(solver.matrix(), sd);
}


void Foam::blockSmootherPrecon::update()
{
    smoother_->update();
}


void Foam::blockSmootherPrecon::precondition
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    w = Zero;
    smoother_->smooth(w, r, nSweeps_);
}


void Foam::blockSmootherPrecon::writeSettings(dictionary& dict) const
{
    dict.add("type", type());
    dict.add("nSweeps", nSweeps_);
}


// ************************************************************************* //
