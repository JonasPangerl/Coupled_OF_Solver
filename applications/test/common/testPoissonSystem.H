/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Description
    Shared set-up of Test-blockMatrix and Test-blockGAMG (spec 6.4):
    a 4x4-block Poisson-like SPD system on the case mesh.

    Scalar operator: M = -laplacian(psi) with fixedValue on all non-constraint
    patches (native fvScalarMatrix; the Dirichlet boundary makes it SPD).
    Block operator: every block is the scalar coefficient times I (4x4):
        diag  = I*D_P            (D = M.D(), including boundary diag)
        upper = I*M.upper_f,  lower = I*M.lower_f
        processor interfaces: I*(-boundaryCoeffs) (native sign -> entry)

    Test vector: analytic, decomposition independent,
        x_k(C) = sin(omega_k . C / L + phi_k),  k = 0..3,
    with L the bounding-box size, so results from 1 and N ranks can be
    compared through global integrals.

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_testPoissonSystem_H
#define coupledFoam_testPoissonSystem_H

#include "fvCFD.H"
#include "blockLduMatrix4.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{
namespace testPoisson
{

//- Scalar field on which the native operator is built
inline tmp<volScalarField> makePsi(const fvMesh& mesh)
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        fixedValueFvPatchScalarField::typeName
    );
}


//- Native scalar operator
inline tmp<fvScalarMatrix> makeScalarMatrix(volScalarField& psi)
{
    return tmp<fvScalarMatrix>::New(-fvm::laplacian(psi));
}


//- Fill the block matrix from the scalar operator
inline void fillBlockMatrix
(
    const fvScalarMatrix& M,
    blockLduMatrix4& B
)
{
    const scalarField D(M.D());
    const scalarField& up = M.upper();
    const scalarField& lo = M.lower();

    B.clear();

    for (label celli = 0; celli < B.nCells(); ++celli)
    {
        blockScalar* d = B.diagBlock(celli);
        for (label k = 0; k < blockDim; ++k)
        {
            d[k*blockDim + k] = narrow(D[celli]);
        }
    }

    for (label facei = 0; facei < B.nFaces(); ++facei)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            B.upper()[facei*blockSize + k*blockDim + k] = narrow(up[facei]);
            B.lower()[facei*blockSize + k*blockDim + k] = narrow(lo[facei]);
        }
    }

    forAll(B.interfaces(), bi)
    {
        const label patchi = B.interfaces()[bi].index();
        const scalarField& bc = M.boundaryCoeffs()[patchi];
        blockScalarList& c = B.interfaceCoeffs(bi);
        forAll(bc, pf)
        {
            for (label k = 0; k < blockDim; ++k)
            {
                c[pf*blockSize + k*blockDim + k] = narrow(-bc[pf]);
            }
        }
    }

    B.markUpdated();
}


//- Analytic test vector component k at point C
inline doubleScalar xExact(const label k, const point& C, const scalar L)
{
    // Wave numbers and phases per component (fixed test data)
    static const doubleScalar omega[4][3] =
    {
        {1.0, 2.0, 0.5},
        {2.5, 0.7, 1.3},
        {0.3, 1.9, 2.2},
        {1.7, 1.1, 0.9}
    };
    static const doubleScalar phi[4] = {0.1, 0.7, 1.3, 2.1};

    const doubleScalar a =
        (omega[k][0]*C.x() + omega[k][1]*C.y() + omega[k][2]*C.z())/L
      + phi[k];
    return std::sin(a);
}


//- Analytic block vector
inline blockScalarList makeX(const fvMesh& mesh)
{
    const scalar L = mag(mesh.bounds().span());
    const volVectorField& C = mesh.C();
    blockScalarList x(blockDim*mesh.nCells());
    forAll(C, celli)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            x[celli*blockDim + k] = narrow(xExact(k, C[celli], L));
        }
    }
    return x;
}

} // End namespace testPoisson
} // End namespace Foam

#endif

// ************************************************************************* //
