/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "MRFCoupling.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFCoupling::MRFCoupling
(
    const fvMesh& mesh,
    const IOMRFZoneList& MRF
)
:
    mesh_(mesh),
    MRF_(MRF),
    coeffs_(),
    active_(MRF.active())
{
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFCoupling::update()
{
    active_ = MRF_.active();

    if (!active_)
    {
        coeffs_.clear();
        return;
    }

    coeffs_.resize_nocopy(mesh_.nCells());
    coeffs_ = Zero;

    const scalarField& V = mesh_.V();

    for (direction j = 0; j < vector::nComponents; ++j)
    {
        vector e(Zero);
        e[j] = 1;

        const volVectorField Ej
        (
            IOobject
            (
                "MRFCoupling::e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedVector(dimVelocity, e)
        );

        // Column j of [Omega]_x is Omega x e_j
        const tmp<volVectorField> tcol = MRF_.DDt(Ej);
        const vectorField& col = tcol().primitiveField();

        forAll(col, celli)
        {
            for (direction i = 0; i < vector::nComponents; ++i)
            {
                coeffs_[celli][i*vector::nComponents + j] = V[celli]*col[celli][i];
            }
        }
    }
}


// ************************************************************************* //
