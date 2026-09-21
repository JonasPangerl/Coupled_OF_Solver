/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGAMGPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockGAMGPrecon, 0);
    addToRunTimeSelectionTable
    (
        blockPreconditioner,
        blockGAMGPrecon,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGAMGPrecon::blockGAMGPrecon
(
    const blockSolver& solver,
    const dictionary& dict
)
:
    blockPreconditioner(solver),
    gamg_(solver.matrix(), dict.subOrEmptyDict("blockGAMG"))
{
    gamg_.writeStats(Info);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockGAMGPrecon::update()
{
    gamg_.update();
}


void Foam::blockGAMGPrecon::precondition
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    gamg_.Vcycle(w, r);
}


void Foam::blockGAMGPrecon::writeSettings(dictionary& dict) const
{
    dict.add("type", type());
    dictionary gd;
    gamg_.writeSettings(gd);
    dict.add("blockGAMG", gd);
}


// ************************************************************************* //
