/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::boundaryCoupling

Description
    Classification of patches for the block assembly (spec 5.5, D-002):

        empty            no faces in the fvMesh (emptyFvPatch), skipped
        processor        implicit block coupling via blockLduInterface
        explicitCoupled  cyclic, cyclicAMI, processorCyclic: neighbour values
                         lagged from the previous iteration (v2 scope)
        physical         every other patch; implicit through the native
                         linearisation U_f = viC . U_P + vbC (and p alike)

    Pressure policy of the Rhie-Chow terms on physical patches (D-013):
    D_f = D_P where p is fixed (fixesValue()), 0 elsewhere.

SourceFiles
    boundaryCoupling.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_boundaryCoupling_H
#define coupledFoam_boundaryCoupling_H

#include "fvPatch.H"
#include "fvPatchFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryCoupling
{

enum class patchKind
{
    empty,
    processor,
    explicitCoupled,
    physical
};

//- Classify a patch
patchKind kind(const fvPatch& patch);

//- Name of a kind (for the log)
const char* kindName(const patchKind k);

//- Per-patch kinds of a mesh
List<patchKind> kinds(const fvMesh& mesh);

//- True if the pressure BC of a physical patch fixes the value
//  (fixedValue, totalPressure, ...), i.e. carries the Rhie-Chow terms
bool fixesPressure(const fvPatchScalarField& pp);

} // End namespace boundaryCoupling
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
