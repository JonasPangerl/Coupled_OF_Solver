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
    volJumpThreshold(scalar(coupledDefaults::volJumpThreshold))
{}


Foam::staticCriteria::settings::settings(const dictionary& staticDict)
:
    settings()
{
    wallStarved = staticDict.getOrDefault<bool>("wallStarved", wallStarved);
    maxWallInternalFaces = staticDict.getOrDefault<label>
    (
        "maxWallInternalFaces",
        maxWallInternalFaces
    );
    procAMI = staticDict.getOrDefault<bool>("procAMI", procAMI);
    volumeJump = staticDict.getOrDefault<bool>("volumeJump", volumeJump);
    volJumpThreshold = staticDict.getOrDefault<scalar>
    (
        "volJumpThreshold",
        volJumpThreshold
    );

    if (maxWallInternalFaces < 0)
    {
        FatalIOErrorInFunction(staticDict)
            << "remediation.static.maxWallInternalFaces must be >= 0, got "
            << maxWallInternalFaces << exit(FatalIOError);
    }
    // Negated comparison also rejects non-finite input
    if (!(volJumpThreshold > 0 && volJumpThreshold < 1))
    {
        FatalIOErrorInFunction(staticDict)
            << "remediation.static.volJumpThreshold must be in (0, 1), got "
            << volJumpThreshold << exit(FatalIOError);
    }
}


void Foam::staticCriteria::settings::write(dictionary& staticDict) const
{
    staticDict.add("wallStarved", wallStarved);
    staticDict.add("maxWallInternalFaces", maxWallInternalFaces);
    staticDict.add("procAMI", procAMI);
    staticDict.add("volumeJump", volumeJump);
    staticDict.add("volJumpThreshold", volJumpThreshold);
}


// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

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
    labelList& bits,
    boolList& isStatic
)
{
    bits.resize_nocopy(mesh.nCells());
    bits = 0;

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

    forAll(bits, celli)
    {
        if (bits[celli])
        {
            if (!isStatic[celli])
            {
                ++c.nNew;
            }
            isStatic[celli] = true;
        }
    }

    // Unconditional on every rank
    reduce(c.nWallStarved, sumOp<label>());
    reduce(c.nProcAMI, sumOp<label>());
    reduce(c.nVolumeJump, sumOp<label>());
    reduce(c.nNew, sumOp<label>());

    return c;
}


void Foam::staticCriteria::report
(
    const settings& s,
    const counts& c,
    const label nTotal
)
{
    // GUARD: nTotal >= 1
    const scalar percentPerCell = scalar(100.0/max(scalar(nTotal), scalar(1)));

    auto item = [&](const bool on, const label n)
    {
        if (on)
        {
            Info<< n << " cells (" << percentPerCell*scalar(n) << " %)";
        }
        else
        {
            Info<< "off";
        }
    };

    Info<< "remediation: static topological criteria (C1): wallStarved"
        << " (<= " << s.maxWallInternalFaces << " internal faces, 3D): ";
    item(s.wallStarved, c.nWallStarved);
    Info<< ", procAMI: ";
    item(s.procAMI, c.nProcAMI);
    Info<< ", volumeJump>" << s.volJumpThreshold << ": ";
    item(s.volumeJump, c.nVolumeJump);
    Info<< "; not marked by the quality criteria: " << c.nNew << " cells"
        << endl;
}


// ************************************************************************* //
