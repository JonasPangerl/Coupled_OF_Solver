/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "zonalSets.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

const Foam::wordList& Foam::zonalSets::builtinNames()
{
    static const wordList names
    ({
        "_wallCells",
        "_procCells",
        "_amiCells",
        "_procAMICells",
        "_wallStarved",
        "_remediationStatic"
    });
    return names;
}


bool Foam::zonalSets::isBuiltin(const wordRe& name)
{
    return name.isLiteral() && builtinNames().found(name);
}


bool Foam::zonalSets::needsStaticSet(const wordRe& name)
{
    return isBuiltin(name) && name == "_remediationStatic";
}


Foam::boolList Foam::zonalSets::mark
(
    const polyMesh& mesh,
    const word& name,
    const staticCriteria::settings& criteria
)
{
    if (name == "_wallCells")
    {
        return staticCriteria::wallCells(mesh);
    }
    if (name == "_procCells")
    {
        return staticCriteria::procCells(mesh);
    }
    if (name == "_amiCells")
    {
        return staticCriteria::amiCells(mesh);
    }
    if (name == "_procAMICells")
    {
        return staticCriteria::procAMI(mesh);
    }
    if (name == "_wallStarved")
    {
        return staticCriteria::wallStarved
        (
            mesh,
            criteria.maxWallInternalFaces
        );
    }

    FatalErrorInFunction
        << "Not a mesh-only built-in zonal set: " << name << nl
        << "Built-in sets: " << flatOutput(builtinNames())
        << exit(FatalError);
    return boolList();
}


// ************************************************************************* //
