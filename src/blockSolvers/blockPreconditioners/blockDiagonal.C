/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockDiagonal.H"
#include "block4Ops.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockDiagonal, 0);
    addToRunTimeSelectionTable(blockPreconditioner, blockDiagonal, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockDiagonal::blockDiagonal
(
    const blockSolver& solver,
    const dictionary& dict
)
:
    blockPreconditioner(solver),
    pivotGuard_
    (
        dict.getOrDefault<doubleScalar>
        (
            "pivotGuard",
            coupledDefaults::pivotGuard
        )
    ),
    rD_(),
    nSingular_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockDiagonal::update()
{
    const blockLduMatrix4& A = matrix();
    const label nCells = A.nCells();

    rD_.resize_nocopy(blockSize*nCells);
    nSingular_ = 0;

    for (label celli = 0; celli < nCells; ++celli)
    {
        nSingular_ += block4Ops::invert
        (
            A.diagBlock(celli),
            rD_.data() + celli*blockSize,
            pivotGuard_
        );
    }
}


void Foam::blockDiagonal::precondition
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    const label nCells = matrix().nCells();

    for (label celli = 0; celli < nCells; ++celli)
    {
        block4Ops::matVec
        (
            rD_.cdata() + celli*blockSize,
            r.cdata() + celli*blockDim,
            w.data() + celli*blockDim
        );
    }
}


void Foam::blockDiagonal::writeSettings(dictionary& dict) const
{
    dict.add("type", type());
    dict.add("pivotGuard", pivotGuard_);
}


// ************************************************************************* //
