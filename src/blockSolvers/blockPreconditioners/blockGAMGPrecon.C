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
    // Amendment B1: the K-cycle is a variable preconditioner
    const word solverType(dict.getOrDefault<word>("solver", word::null));
    if
    (
        gamg_.cycleType() == blockGAMG::cycleKind::K
     && solverType != "blockFGMRES"
    )
    {
        FatalIOErrorInFunction(dict)
            << "blockGAMG cycleType K requires solver blockFGMRES (found "
            << solverType << "): the K-cycle is a variable preconditioner"
            << " (amendment B1)" << exit(FatalIOError);
    }

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
    gamg_.apply(w, r);
}


void Foam::blockGAMGPrecon::writeSettings(dictionary& dict) const
{
    dict.add("type", type());
    dictionary gd;
    gamg_.writeSettings(gd);
    dict.add("blockGAMG", gd);
}


// ************************************************************************* //
