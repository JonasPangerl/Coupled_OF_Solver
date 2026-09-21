/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockPreconditioner, 0);
    defineRunTimeSelectionTable(blockPreconditioner, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blockPreconditioner> Foam::blockPreconditioner::New
(
    const blockSolver& solver,
    const dictionary& dict
)
{
    const word preconType(dict.get<word>("preconditioner"));

    auto* ctorPtr = dictionaryConstructorTable(preconType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "blockPreconditioner",
            preconType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<blockPreconditioner>(ctorPtr(solver, dict));
}


// ************************************************************************* //
