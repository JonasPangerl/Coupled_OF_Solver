/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "boundaryCoupling.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

Foam::boundaryCoupling::patchKind Foam::boundaryCoupling::kind
(
    const fvPatch& patch
)
{
    if (patch.size() == 0 && !patch.coupled())
    {
        return patchKind::empty;
    }
    if (isA<processorCyclicFvPatch>(patch))
    {
        // Transformation across the interface not handled by the block
        // exchange: lagged like cyclic. blockLduInterface::isBlockCoupled
        // excludes these interfaces accordingly (no block interface, no
        // exchange in Amul).
        return patchKind::explicitCoupled;
    }
    if (isA<processorFvPatch>(patch))
    {
        return patchKind::processor;
    }
    if (patch.coupled())
    {
        return patchKind::explicitCoupled;
    }
    return patchKind::physical;
}


const char* Foam::boundaryCoupling::kindName(const patchKind k)
{
    switch (k)
    {
        case patchKind::empty: return "empty";
        case patchKind::processor: return "processor";
        case patchKind::explicitCoupled: return "explicitCoupled";
        case patchKind::physical: return "physical";
    }
    return "unknown";
}


Foam::List<Foam::boundaryCoupling::patchKind>
Foam::boundaryCoupling::kinds(const fvMesh& mesh)
{
    List<patchKind> k(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        k[patchi] = kind(mesh.boundary()[patchi]);
    }
    return k;
}


bool Foam::boundaryCoupling::fixesPressure(const fvPatchScalarField& pp)
{
    return pp.fixesValue();
}


// ************************************************************************* //
