/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockSmoother.H"
#include "coupledDefaults.H"
#include "blockScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockSmoother, 0);
    defineRunTimeSelectionTable(blockSmoother, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockSmoother::blockSmoother
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    matrix_(matrix),
    pivotGuard_
    (
        dict.getOrDefault<doubleScalar>
        (
            "pivotGuard",
            coupledDefaults::pivotGuard
        )
    ),
    nSingular_(0),
    relax_
    (
        narrow
        (
            dict.getOrDefault<doubleScalar>
            (
                "smootherRelaxation",
                coupledDefaults::smootherRelaxation
            )
        )
    ),
    damped_(relax_ != blockScalar(1))
{
    if (!(relax_ > 0) || relax_ > 1)
    {
        FatalIOErrorInFunction(dict)
            << "smootherRelaxation must be in (0, 1], got " << relax_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blockSmoother> Foam::blockSmoother::New
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
{
    const word smootherType
    (
        dict.getOrDefault<word>("smoother", word(coupledDefaults::smoother))
    );

    auto* ctorPtr = dictionaryConstructorTable(smootherType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "blockSmoother",
            smootherType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<blockSmoother>(ctorPtr(matrix, dict));
}


// ************************************************************************* //
