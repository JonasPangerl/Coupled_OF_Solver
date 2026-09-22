/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::staticCriteria

Description
    Topological static-set criteria of amendment C1 (spec 8.1 extended).
    Pure functions of the mesh and the thresholds; they return per-cell
    markers and bit masks, and the caller (remediation::buildStatic) ORs
    them into the static set. Evaluated once at start-up.

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

    Bits in remediationFlag (flags are OR-ed): 1 static (any static
    criterion), 2 dynamic, 4 wallStarved, 8 procAMI, 16 volumeJump.

    Dictionary (coupled.remediation.static), defaults in coupledDefaults.H:
    \verbatim
    wallStarved             yes;    // switch (coupledFoam addition)
    maxWallInternalFaces    2;
    procAMI                 yes;
    volumeJump              yes;    // switch (coupledFoam addition)
    volJumpThreshold        0.98;   // i.e. volume ratio > 50
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

//- remediationFlag bits (spec 8.1/12, amendment C1)
enum flagBits : label
{
    bitStatic = 1,
    bitDynamic = 2,
    bitWallStarved = 4,
    bitProcAMI = 8,
    bitVolumeJump = 16
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

    //- Read from coupled.remediation.static (checks the input)
    explicit settings(const dictionary& staticDict);

    //- Effective settings into the static sub-dictionary
    void write(dictionary& staticDict) const;
};


//- Global numbers of cells marked by each criterion
struct counts
{
    label nWallStarved = 0;
    label nProcAMI = 0;
    label nVolumeJump = 0;
    //- Cells marked by a topological criterion only (not by 8.1)
    label nNew = 0;
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

//- Evaluate the enabled criteria: fills bits (per cell, 4/8/16), ORs the
//  marked cells into isStatic and returns the global counts. isStatic must
//  already hold the 8.1 quality criteria (for counts.nNew). Collective.
counts apply
(
    const polyMesh& mesh,
    const settings& s,
    labelList& bits,
    boolList& isStatic
);

//- One-line summary for the once-per-run block (12.2)
void report(const settings& s, const counts& c, const label nTotal);


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
