/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "nonOrthCorrection.H"
#include "snGradScheme.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::nonOrthCorrection::limitedCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const scalar lambda
) const
{
    // Native "limited corrected <lambda>" scheme (limitedSnGrad)
    IStringStream schemeData("limited corrected " + Foam::name(lambda));
    tmp<fv::snGradScheme<Type>> tscheme
    (
        fv::snGradScheme<Type>::New(mesh_, schemeData)
    );
    return tscheme().correction(vf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::nonOrthCorrection::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfFieldType;

    if (orthogonal_)
    {
        return tmp<surfFieldType>::New
        (
            IOobject
            (
                "nonOrthCorrection(" + vf.name() + ')',
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensioned<Type>(vf.dimensions()/dimLength, Zero)
        );
    }

    tmp<surfFieldType> tc = limitedCorrection(vf, limiter_);

    if (!anyStatic_ || limiterStatic_ == limiter_)
    {
        return tc;
    }

    const tmp<surfFieldType> tcs = limitedCorrection(vf, limiterStatic_);

    surfFieldType& c = tc.ref();
    const surfFieldType& cs = tcs();

    Field<Type>& ci = c.primitiveFieldRef();
    const Field<Type>& csi = cs.primitiveField();
    forAll(staticFace_, facei)
    {
        if (staticFace_[facei])
        {
            ci[facei] = csi[facei];
        }
    }

    auto& cbf = c.boundaryFieldRef();
    forAll(cbf, patchi)
    {
        const boolList& sp = staticPatchFace_[patchi];
        forAll(sp, pf)
        {
            if (sp[pf])
            {
                cbf[patchi][pf] = cs.boundaryField()[patchi][pf];
            }
        }
    }

    return tc;
}


// ************************************************************************* //
