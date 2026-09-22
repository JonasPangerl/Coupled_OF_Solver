/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::zonalSets

Description
    Built-in named cell sets of the zonal factors (amendment C5, spec 8.3
    extended), usable in place of a cellZone name in
    coupled.remediation.zonal.zones:
    \verbatim
    zones
    (
        { cellZone _wallCells;          cflFactor 0.7; }
        { cellZone _remediationStatic;  beta 0.5; }
    );
    \endverbatim

    Names (literal, not regular expressions):
      _wallCells           cells with a wall face
      _procCells           cells with a processor face (decomposition
                           dependent)
      _amiCells            cells with a cyclicAMI face
      _procAMICells        cells with a processor and a cyclicAMI face
                           (= the C1 procAMI criterion)
      _wallStarved         the C1 wallStarved criterion with the
                           remediation.static.maxWallInternalFaces setting
                           (evaluated even if the static switch is off)
      _remediationStatic   the complete static set (8.1 + C1); applied once
                           the static set is built (remediation::buildStatic)

    A literal name that equals a built-in name always means the built-in
    set; regular expressions only match real cellZones. The patchDistance
    patches were already native wordRes (regular expressions and patch
    groups, B6). A cell affected by several entries gets the product of the
    factors (unchanged rule of D-027).

SourceFiles
    zonalSets.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_zonalSets_H
#define coupledFoam_zonalSets_H

#include "polyMesh.H"
#include "wordRe.H"
#include "boolList.H"
#include "staticCriteria.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zonalSets
{

//- A zonal entry whose cells are only known after the static set is built
struct deferredEntry
{
    word name;
    scalar cflFactor;
    scalar beta;
};

//- All built-in names
const wordList& builtinNames();

//- Is the (literal) name a built-in set?
bool isBuiltin(const wordRe& name);

//- Does the built-in set need the static set (_remediationStatic)?
bool needsStaticSet(const wordRe& name);

//- Cells of a mesh-only built-in set (not _remediationStatic). Local.
boolList mark
(
    const polyMesh& mesh,
    const word& name,
    const staticCriteria::settings& criteria
);


} // End namespace zonalSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
