/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "staticCriteria.H"
#include "coupledDefaults.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "wallPolyPatch.H"
#include "syncTools.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- Patch with an implicit neighbour cell (processor or cyclic non-AMI)
bool implicitPatch(const Foam::polyPatch& pp)
{
    return
        Foam::isA<Foam::processorPolyPatch>(pp)
     || Foam::isA<Foam::cyclicPolyPatch>(pp);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * settings  * * * * * * * * * * * * * * * //

Foam::staticCriteria::settings::settings()
:
    wallStarved(coupledDefaults::wallStarvedEnabled),
    maxWallInternalFaces(coupledDefaults::maxWallInternalFaces),
    procAMI(coupledDefaults::procAMIEnabled),
    volumeJump(coupledDefaults::volumeJumpEnabled),
    volJumpThreshold(coupledDefaults::volJumpThreshold)
{}


void Foam::staticCriteria::settings::check(const dictionary& errDict) const
{
    if (maxWallInternalFaces < 0)
    {
        FatalIOErrorInFunction(errDict)
            << "remediation.wall.maxWallInternalFaces must be >= 0, got "
            << maxWallInternalFaces << exit(FatalIOError);
    }
    // Negated comparison also rejects non-finite input
    if (!(volJumpThreshold > 0 && volJumpThreshold < 1))
    {
        FatalIOErrorInFunction(errDict)
            << "remediation.badMesh.volJumpThreshold must be in (0, 1), got "
            << volJumpThreshold << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

Foam::scalar Foam::staticCriteria::limitedNonOrthCoeff
(
    const dictionary& coupledDict
)
{
    const dictionary& r = coupledDict.subOrEmptyDict("remediation");
    const scalar legacy =
        r.subOrEmptyDict("static").getOrDefault<scalar>
        (
            "nonOrthLimiter",
            coupledDefaults::limitedNonOrthCoeff
        );
    const scalar c = r.getOrDefault<scalar>("limitedNonOrthCoeff", legacy);
    if
    (
        r.found("limitedNonOrthCoeff")
     && r.subOrEmptyDict("static").found("nonOrthLimiter")
     && !(c == legacy)
    )
    {
        FatalIOErrorInFunction(r)
            << "remediation.limitedNonOrthCoeff and the deprecated"
            << " remediation.static.nonOrthLimiter are both given, with"
            << " different values (D-066)" << exit(FatalIOError);
    }
    // Negated comparison also rejects non-finite input
    if (!(c >= 0 && c <= 1))
    {
        FatalIOErrorInFunction(r)
            << "remediation.limitedNonOrthCoeff must be in [0, 1], got " << c
            << exit(FatalIOError);
    }
    return c;
}


Foam::string Foam::staticCriteria::limitedGradScheme
(
    const dictionary& coupledDict
)
{
    return coupledDict.subOrEmptyDict("remediation").getOrDefault<string>
    (
        "limitedGradScheme",
        string(coupledDefaults::limitedGradScheme)
    );
}


bool Foam::staticCriteria::hasImplicitNeighbour
(
    const polyMesh& mesh,
    const label facei
)
{
    if (facei < mesh.nInternalFaces())
    {
        return true;
    }
    const label patchi = mesh.boundaryMesh().whichPatch(facei);
    return patchi >= 0 && implicitPatch(mesh.boundaryMesh()[patchi]);
}


Foam::boolList Foam::staticCriteria::wallStarved
(
    const polyMesh& mesh,
    const label maxInternal
)
{
    boolList mark(mesh.nCells(), false);

    // 3D cases only (C1)
    if (mesh.nSolutionD() != vector::nComponents)
    {
        return mark;
    }

    labelList nInternal(mesh.nCells(), Zero);
    boolList hasWall(mesh.nCells(), false);

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        ++nInternal[own[facei]];
        ++nInternal[nei[facei]];
    }

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        const bool implicit = implicitPatch(pp);
        const bool wall = isA<wallPolyPatch>(pp);
        for (const label celli : pp.faceCells())
        {
            if (implicit) ++nInternal[celli];
            if (wall) hasWall[celli] = true;
        }
    }

    forAll(mark, celli)
    {
        mark[celli] = hasWall[celli] && nInternal[celli] <= maxInternal;
    }
    return mark;
}


Foam::boolList Foam::staticCriteria::wallCells(const polyMesh& mesh)
{
    return cellsOnPatches
    (
        mesh,
        [](const polyPatch& pp) { return isA<wallPolyPatch>(pp); }
    );
}


Foam::boolList Foam::staticCriteria::procCells(const polyMesh& mesh)
{
    return cellsOnPatches
    (
        mesh,
        [](const polyPatch& pp) { return isA<processorPolyPatch>(pp); }
    );
}


Foam::boolList Foam::staticCriteria::amiCells(const polyMesh& mesh)
{
    return cellsOnPatches
    (
        mesh,
        [](const polyPatch& pp) { return isA<cyclicAMIPolyPatch>(pp); }
    );
}


Foam::boolList Foam::staticCriteria::procAMI(const polyMesh& mesh)
{
    boolList mark(procCells(mesh));
    const boolList ami(amiCells(mesh));
    forAll(mark, celli)
    {
        mark[celli] = mark[celli] && ami[celli];
    }
    return mark;
}


Foam::boolList Foam::staticCriteria::volumeJump
(
    const polyMesh& mesh,
    const scalar threshold
)
{
    boolList mark(mesh.nCells(), false);

    const scalarField& V = mesh.cellVolumes();

    // |a - b| > threshold*max(a, b), written without a division
    auto jump = [threshold](const scalar a, const scalar b)
    {
        const scalar m = max(a, b);
        return m > 0 && mag(a - b) > threshold*m;
    };

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        if (jump(V[own[facei]], V[nei[facei]]))
        {
            mark[own[facei]] = true;
            mark[nei[facei]] = true;
        }
    }

    // Neighbour volumes across processor and cyclic faces (collective)
    scalarList Vnbr;
    syncTools::swapBoundaryCellList(mesh, V, Vnbr);

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (!implicitPatch(pp))
        {
            continue;
        }
        const labelUList& fc = pp.faceCells();
        const label offset = pp.offset();
        forAll(fc, i)
        {
            if (jump(V[fc[i]], Vnbr[offset + i]))
            {
                mark[fc[i]] = true;
            }
        }
    }

    return mark;
}


Foam::staticCriteria::counts Foam::staticCriteria::apply
(
    const polyMesh& mesh,
    const settings& s,
    labelList& bits
)
{
    counts c;

    auto add = [&bits](const boolList& mark, const label bit, label& n)
    {
        forAll(mark, celli)
        {
            if (mark[celli])
            {
                bits[celli] |= bit;
                ++n;
            }
        }
    };

    // The switches are identical on all ranks: the collective inside
    // volumeJump() is either called everywhere or nowhere
    if (s.wallStarved)
    {
        add
        (
            wallStarved(mesh, s.maxWallInternalFaces),
            bitWallStarved,
            c.nWallStarved
        );
    }
    if (s.procAMI)
    {
        add(procAMI(mesh), bitProcAMI, c.nProcAMI);
    }
    if (s.volumeJump)
    {
        add
        (
            volumeJump(mesh, s.volJumpThreshold),
            bitVolumeJump,
            c.nVolumeJump
        );
    }

    // Unconditional on every rank
    reduce(c.nWallStarved, sumOp<label>());
    reduce(c.nProcAMI, sumOp<label>());
    reduce(c.nVolumeJump, sumOp<label>());

    return c;
}


// ************************************************************************* //
