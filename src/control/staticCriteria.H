/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::staticCriteria

Description
    Topological static-set criteria of amendment C1 (spec 8.1 extended).
    Pure functions of the mesh and the thresholds; they return per-cell
    markers and bit masks, and the caller (remediation::buildStatic) sorts
    them into the pre-selected cell categories (D-066: wallStarved -> wall,
    procAMI -> processor, volumeJump -> badMesh). Evaluated once at
    start-up.

    - wallStarved: a cell with at least one wall face and at most
      maxWallInternalFaces "internal" faces (3D cases only, i.e.
      mesh.nSolutionD() == 3). Internal means a face with an implicit
      neighbour cell: an internal face, a processor face (an internal face
      of the undecomposed mesh; this keeps the criterion independent of the
      decomposition) or a cyclic (non-AMI) face. cyclicAMI faces are not
      counted (their coupling is lagged, spec 5.5).
    - procAMI: a cell with both a processor face and a cyclicAMI face. By
      construction this set depends on the decomposition (it grows with the
      number of ranks).
    - volumeJump: |V_P - V_N|/max(V_P, V_N) > volJumpThreshold across an
      internal, processor or cyclic (non-AMI) face; both cells are marked.
      Complements volRatioThreshold (8.1).

    Bits in remediationFlag (flags are OR-ed; enum flagBits below, D-066):
      1 static (any pre-selected category), 2 dynamic,
      4 wallStarved, 8 procAMI, 16 volumeJump (criteria, C1),
      32 meshQuality, 64 badMesh, 128 processor, 256 wall (categories),
      512 closedness, 1024 severe non-orthogonality/skewness,
      2048 processor-patch layer, 4096 wall-patch layer (criteria).

    Keywords (coupled.remediation.<category>, D-066; the pre-D-066
    remediation.static entries map onto them), defaults in
    coupledDefaults.H:
    \verbatim
    wall      { wallStarved no; maxWallInternalFaces 2; }
    processor { procAMI yes; }
    badMesh   { volumeJump yes; volJumpThreshold 0.98; }  // ratio > 50
    \endverbatim

    Collectives: volumeJump() swaps the neighbour-cell volumes over the
    coupled boundaries (syncTools, collective) and counts() reduces; both
    must be called on every rank. wallStarved() and procAMI() are local.

SourceFiles
    staticCriteria.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_staticCriteria_H
#define coupledFoam_staticCriteria_H

#include "polyMesh.H"
#include "boolList.H"
#include "labelList.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace staticCriteria
{

//- remediationFlag bits (spec 8.1/12, amendment C1, categories D-066)
enum flagBits : label
{
    bitStatic = 1,               // any pre-selected category
    bitDynamic = 2,
    // C1 criteria
    bitWallStarved = 4,
    bitProcAMI = 8,
    bitVolumeJump = 16,
    // Categories (D-066)
    bitMeshQuality = 32,
    bitBadMesh = 64,
    bitProcessor = 128,
    bitWall = 256,
    // Further criteria (D-066)
    bitClosedness = 512,
    bitSevereQuality = 1024,
    bitProcLayer = 2048,
    bitWallLayer = 4096
};


//- Settings of the topological criteria
struct settings
{
    bool wallStarved;
    label maxWallInternalFaces;
    bool procAMI;
    bool volumeJump;
    scalar volJumpThreshold;

    //- Defaults from coupledDefaults.H
    settings();

    //- Check the values; errors refer to the given dictionary
    void check(const dictionary& errDict) const;
};


//- Global numbers of cells marked by each criterion
struct counts
{
    label nWallStarved = 0;
    label nProcAMI = 0;
    label nVolumeJump = 0;
};


//- Does the face have an implicit neighbour cell? (internal, processor or
//  cyclic non-AMI face)
bool hasImplicitNeighbour(const polyMesh& mesh, const label facei);

//- Cells with a wall face and <= maxInternal implicit-neighbour faces
//  (empty marker if mesh.nSolutionD() != 3). Local.
boolList wallStarved(const polyMesh& mesh, const label maxInternal);

//- Cells with at least one face on a patch for which pred(patch) is true.
//  Local.
template<class Predicate>
boolList cellsOnPatches(const polyMesh& mesh, const Predicate& pred);

//- Cells with a wall face. Local.
boolList wallCells(const polyMesh& mesh);

//- Cells with a processor face. Local.
boolList procCells(const polyMesh& mesh);

//- Cells with a cyclicAMI face. Local.
boolList amiCells(const polyMesh& mesh);

//- Cells with both a processor face and a cyclicAMI face. Local.
boolList procAMI(const polyMesh& mesh);

//- Cells across a face with |V_P - V_N|/max(V_P, V_N) > threshold.
//  Collective (neighbour volumes over coupled boundaries).
boolList volumeJump(const polyMesh& mesh, const scalar threshold);

//- Non-orthogonal limiter coefficient of the cells with nonOrthLimiter
//  (coupled.remediation.limitedNonOrthCoeff; legacy
//  coupled.remediation.static.nonOrthLimiter; D-066). Checked: [0, 1].
scalar limitedNonOrthCoeff(const dictionary& coupledDict);

//- Gradient scheme of p in the Rhie-Chow terms of the cells with
//  gradLimiter (coupled.remediation.limitedGradScheme, D-018/D-066)
string limitedGradScheme(const dictionary& coupledDict);

//- Evaluate the enabled criteria: ORs bits 4/8/16 into the per-cell bits
//  (sized nCells) and returns the global counts. Collective.
counts apply
(
    const polyMesh& mesh,
    const settings& s,
    labelList& bits
);


} // End namespace staticCriteria
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Predicate>
Foam::boolList Foam::staticCriteria::cellsOnPatches
(
    const polyMesh& mesh,
    const Predicate& pred
)
{
    boolList mark(mesh.nCells(), false);
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (pred(pp))
        {
            for (const label celli : pp.faceCells())
            {
                mark[celli] = true;
            }
        }
    }
    return mark;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
