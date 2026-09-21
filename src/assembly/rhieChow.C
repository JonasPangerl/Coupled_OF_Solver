/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "rhieChow.H"
#include "nonOrthCorrection.H"
#include "fvcGrad.H"
#include "calculatedFvPatchFields.H"
#include "calculatedFvsPatchFields.H"
#include "gradScheme.H"
#include "IStringStream.H"
#include "PstreamReduceOps.H"

using namespace Foam::boundaryCoupling;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhieChow::rhieChow(const fvMesh& mesh)
:
    mesh_(mesh),
    kinds_(boundaryCoupling::kinds(mesh)),
    D_
    (
        IOobject
        (
            "rhieChowD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimTime, Zero),
        calculatedFvPatchScalarField::typeName
    ),
    Df_
    (
        IOobject
        (
            "rhieChowDf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimTime, Zero),
        calculatedFvsPatchScalarField::typeName
    ),
    gradp_
    (
        IOobject
        (
            "rhieChowGradp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero),
        calculatedFvPatchVectorField::typeName
    ),
    q_
    (
        IOobject
        (
            "rhieChowQ",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Zero),
        calculatedFvsPatchScalarField::typeName
    ),
    valid_(false),
    isStatic_(mesh.nCells(), false),
    anyStatic_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhieChow::updateD
(
    const scalarField& abar,
    const volScalarField& p
)
{
    const scalarField& V = mesh_.V();
    scalarField& Di = D_.primitiveFieldRef();

    forAll(Di, celli)
    {
        // GUARD: abar > 0 for a PTC-augmented momentum diagonal; guard anyway
        Di[celli] = V[celli]/max(abar[celli], VSMALL);
    }
    D_.correctBoundaryConditions();

    const surfaceScalarField& w = mesh_.weights();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    scalarField& Dfi = Df_.primitiveFieldRef();
    forAll(Dfi, facei)
    {
        const scalar wf = w[facei];
        Dfi[facei] = wf*Di[own[facei]] + (1 - wf)*Di[nei[facei]];
    }

    auto& Dfb = Df_.boundaryFieldRef();
    forAll(Dfb, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        scalarField& Dfp = Dfb[patchi];
        const labelUList& fc = patch.faceCells();

        switch (kinds_[patchi])
        {
            case patchKind::empty:
                break;

            case patchKind::processor:
            case patchKind::explicitCoupled:
            {
                const scalarField& wp = w.boundaryField()[patchi];
                const scalarField Dn
                (
                    D_.boundaryField()[patchi].patchNeighbourField()
                );
                forAll(Dfp, pf)
                {
                    Dfp[pf] = wp[pf]*Di[fc[pf]] + (1 - wp[pf])*Dn[pf];
                }
                break;
            }

            case patchKind::physical:
            {
                const bool fixedP =
                    fixesPressure(p.boundaryField()[patchi]);
                forAll(Dfp, pf)
                {
                    Dfp[pf] = fixedP ? Di[fc[pf]] : 0;
                }
                break;
            }
        }
    }

    valid_ = true;
}


void Foam::rhieChow::setStaticCells(const boolList& isStatic)
{
    isStatic_ = isStatic;
    bool any = false;
    forAll(isStatic_, celli)
    {
        any = any || isStatic_[celli];
    }
    anyStatic_ = returnReduceOr(any);
}


Foam::scalar Foam::rhieChow::Dref() const
{
    return gAverage(D_.primitiveField());
}


Foam::scalar Foam::rhieChow::g(const label facei) const
{
    return
        Df_[facei]*mesh_.magSf()[facei]
       *mesh_.nonOrthDeltaCoeffs()[facei];
}


Foam::scalar Foam::rhieChow::g(const label patchi, const label pf) const
{
    return
        Df_.boundaryField()[patchi][pf]
       *mesh_.magSf().boundaryField()[patchi][pf]
       *mesh_.nonOrthDeltaCoeffs().boundaryField()[patchi][pf];
}


void Foam::rhieChow::updateExplicit
(
    const volScalarField& p,
    const nonOrthCorrection& noc
)
{
    gradp_ = fvc::grad(p);

    if (anyStatic_)
    {
        // Static remediation cells: limited gradient (spec 8.1, D-018)
        IStringStream schemeData("cellLimited Gauss linear 1");
        tmp<fv::gradScheme<scalar>> tscheme =
            fv::gradScheme<scalar>::New(mesh_, schemeData);
        const tmp<volVectorField> tgl = tscheme().grad(p, "grad(p)Limited");
        const vectorField& gl = tgl().primitiveField();
        vectorField& g = gradp_.primitiveFieldRef();
        forAll(g, celli)
        {
            if (isStatic_[celli])
            {
                g[celli] = gl[celli];
            }
        }
    }

    gradp_.correctBoundaryConditions();

    // Limited non-orthogonal correction from the same gradient, so the
    // cellLimited gradient of static cells (D-018) enters it as well
    const tmp<surfaceScalarField> tcorr = noc.correctionFromGrad(p, gradp_);
    const surfaceScalarField& corr = tcorr();

    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& gp = gradp_.primitiveField();

    scalarField& qi = q_.primitiveFieldRef();
    forAll(qi, facei)
    {
        const scalar wf = w[facei];
        const vector gf = wf*gp[own[facei]] + (1 - wf)*gp[nei[facei]];
        qi[facei] =
            Df_[facei]*((gf & Sf[facei]) - magSf[facei]*corr[facei]);
    }

    auto& qb = q_.boundaryFieldRef();
    forAll(qb, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        scalarField& qp = qb[patchi];
        const labelUList& fc = patch.faceCells();
        const vectorField& Sp = Sf.boundaryField()[patchi];
        const scalarField& Dfp = Df_.boundaryField()[patchi];

        switch (kinds_[patchi])
        {
            case patchKind::empty:
                break;

            case patchKind::processor:
            case patchKind::explicitCoupled:
            {
                const scalarField& wp = w.boundaryField()[patchi];
                const scalarField& mp = magSf.boundaryField()[patchi];
                const scalarField& cp = corr.boundaryField()[patchi];
                const vectorField gn
                (
                    gradp_.boundaryField()[patchi].patchNeighbourField()
                );
                forAll(qp, pf)
                {
                    const vector gf =
                        wp[pf]*gp[fc[pf]] + (1 - wp[pf])*gn[pf];
                    qp[pf] = Dfp[pf]*((gf & Sp[pf]) - mp[pf]*cp[pf]);
                }
                break;
            }

            case patchKind::physical:
            {
                // The implicit p-p term of the face is g (1 - viC_p)
                // (p_b = viC p_P + vbC); the explicit gradient term carries
                // the same weight, so that on a mixed p condition (all mixed
                // types report fixesValue()) a face in zero-gradient mode
                // (valueFraction 0, viC 1) gets no Rhie-Chow flux at all.
                // fixedValue: viC = 0, weight 1; zeroGradient: D_f = 0.
                const scalarField viCp
                (
                    p.boundaryField()[patchi].valueInternalCoeffs
                    (
                        w.boundaryField()[patchi]
                    )
                );
                forAll(qp, pf)
                {
                    qp[pf] =
                        (1 - viCp[pf])*Dfp[pf]*(gp[fc[pf]] & Sp[pf]);
                }
                break;
            }
        }
    }
}


void Foam::rhieChow::updateFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const nonOrthCorrection& noc,
    const bool recomputeExplicit
)
{
    if (recomputeExplicit)
    {
        updateExplicit(p, noc);
    }

    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();

    scalarField& phii = phi.primitiveFieldRef();
    forAll(phii, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wf = w[facei];
        phii[facei] =
            ((wf*Ui[P] + (1 - wf)*Ui[N]) & Sf[facei])
          - g(facei)*(pi[N] - pi[P])
          + q_[facei];
    }

    auto& phib = phi.boundaryFieldRef();
    forAll(phib, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        scalarField& phip = phib[patchi];
        const labelUList& fc = patch.faceCells();
        const vectorField& Sp = Sf.boundaryField()[patchi];
        const scalarField& qp = q_.boundaryField()[patchi];

        switch (kinds_[patchi])
        {
            case patchKind::empty:
                break;

            case patchKind::processor:
            case patchKind::explicitCoupled:
            {
                const scalarField& wp = w.boundaryField()[patchi];
                const vectorField Un
                (
                    U.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField pn
                (
                    p.boundaryField()[patchi].patchNeighbourField()
                );
                forAll(phip, pf)
                {
                    const label P = fc[pf];
                    phip[pf] =
                        ((wp[pf]*Ui[P] + (1 - wp[pf])*Un[pf]) & Sp[pf])
                      - g(patchi, pf)*(pn[pf] - pi[P])
                      + qp[pf];
                }
                break;
            }

            case patchKind::physical:
            {
                const vectorField& Ub = U.boundaryField()[patchi];
                const scalarField& pb = p.boundaryField()[patchi];
                forAll(phip, pf)
                {
                    const label P = fc[pf];
                    phip[pf] =
                        (Ub[pf] & Sp[pf])
                      - g(patchi, pf)*(pb[pf] - pi[P])
                      + qp[pf];
                }
                break;
            }
        }
    }
}


void Foam::rhieChow::writeFields() const
{
    surfaceScalarField Df
    (
        IOobject
        (
            "rhieChowDf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        Df_
    );
    Df.write();
}


// ************************************************************************* //
