/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "sentinel.H"
#include "coupledDefaults.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "globalIndex.H"
#include "fileOperation.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sentinel::sentinel(const fvMesh& mesh, const dictionary& coupledDict)
:
    mesh_(mesh),
    maxRollbacks_(coupledDefaults::maxRollbacks),
    UFactor_(coupledDefaults::sentinelUFactor),
    pFactor_(coupledDefaults::sentinelPFactor),
    cflFactor_(coupledDefaults::sentinelCflFactor),
    maxReport_(coupledDefaults::sentinelMaxReport),
    consecutive_(0),
    nRollbacks_(0),
    stored_(false)
{
    const dictionary& d = coupledDict.subOrEmptyDict("sentinel");
    maxRollbacks_ = d.getOrDefault<label>("maxRollbacks", maxRollbacks_);
    UFactor_ = d.getOrDefault<scalar>("UFactor", UFactor_);
    pFactor_ = d.getOrDefault<scalar>("pFactor", pFactor_);
    cflFactor_ = d.getOrDefault<scalar>("cflFactor", cflFactor_);
    maxReport_ = d.getOrDefault<label>("maxReport", maxReport_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sentinel::store
(
    const volVectorField& U,
    const volScalarField& p,
    const surfaceScalarField& phi,
    const volScalarField* k,
    const volScalarField* omega,
    const volScalarField* nut
)
{
    U0_ = U.primitiveField();
    p0_ = p.primitiveField();
    phi0_ = phi.primitiveField();
    phi0Boundary_.resize(phi.boundaryField().size());
    forAll(phi.boundaryField(), patchi)
    {
        phi0Boundary_[patchi] = phi.boundaryField()[patchi];
    }
    if (k) k0_ = k->primitiveField();
    if (omega) omega0_ = omega->primitiveField();
    if (nut) nut0_ = nut->primitiveField();
    stored_ = true;
}


Foam::sentinel::checkResult Foam::sentinel::check
(
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField* k,
    const volScalarField* omega,
    const scalar Uref,
    const scalar pref
) const
{
    checkResult res;

    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();

    scalar maxU = 0;
    scalar minP = GREAT, maxP = -GREAT;
    scalar minK = GREAT, maxK = -GREAT;
    scalar minW = GREAT, maxW = -GREAT;
    label nNonFinite = 0;

    const scalar Ulim = UFactor_*Uref;
    const scalar plim = pFactor_*pref;

    DynamicList<label> bad;

    forAll(Ui, celli)
    {
        bool cellBad = false;

        const vector& u = Ui[celli];
        if
        (
            !std::isfinite(u.x()) || !std::isfinite(u.y())
         || !std::isfinite(u.z()) || !std::isfinite(pi[celli])
        )
        {
            ++nNonFinite;
            cellBad = true;
        }
        else
        {
            const scalar mu = mag(u);
            maxU = max(maxU, mu);
            minP = min(minP, pi[celli]);
            maxP = max(maxP, pi[celli]);
            cellBad = (mu > Ulim) || (mag(pi[celli]) > plim);
        }

        if (k)
        {
            const scalar kv = (*k)[celli];
            if (!std::isfinite(kv)) { ++nNonFinite; cellBad = true; }
            else { minK = min(minK, kv); maxK = max(maxK, kv); }
        }
        if (omega)
        {
            const scalar wv = (*omega)[celli];
            if (!std::isfinite(wv)) { ++nNonFinite; cellBad = true; }
            else { minW = min(minW, wv); maxW = max(maxW, wv); }
        }

        if (cellBad)
        {
            bad.append(celli);
        }
    }

    reduce(maxU, maxOp<scalar>());
    reduce(minP, minOp<scalar>());
    reduce(maxP, maxOp<scalar>());
    reduce(minK, minOp<scalar>());
    reduce(maxK, maxOp<scalar>());
    reduce(minW, minOp<scalar>());
    reduce(maxW, maxOp<scalar>());
    reduce(nNonFinite, sumOp<label>());

    res.nNonFinite = nNonFinite;
    res.maxMagU = maxU;
    res.minP = minP;
    res.maxP = maxP;
    res.minK = minK;
    res.maxK = maxK;
    res.minOmega = minW;
    res.maxOmega = maxW;
    res.ok =
        nNonFinite == 0
     && maxU <= Ulim
     && max(mag(minP), mag(maxP)) <= plim;
    res.offending = labelList(std::move(bad));

    return res;
}


void Foam::sentinel::restore
(
    volVectorField& U,
    volScalarField& p,
    surfaceScalarField& phi,
    volScalarField* k,
    volScalarField* omega,
    volScalarField* nut
)
{
    if (!stored_)
    {
        FatalErrorInFunction
            << "No stored fields to roll back to" << exit(FatalError);
    }

    U.primitiveFieldRef() = U0_;
    p.primitiveFieldRef() = p0_;
    phi.primitiveFieldRef() = phi0_;
    forAll(phi.boundaryField(), patchi)
    {
        phi.boundaryFieldRef()[patchi] == phi0Boundary_[patchi];
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    if (k && k0_.size())
    {
        k->primitiveFieldRef() = k0_;
        k->correctBoundaryConditions();
    }
    if (omega && omega0_.size())
    {
        omega->primitiveFieldRef() = omega0_;
        omega->correctBoundaryConditions();
    }
    if (nut && nut0_.size())
    {
        nut->primitiveFieldRef() = nut0_;
        nut->correctBoundaryConditions();
    }

    ++consecutive_;
    ++nRollbacks_;
}


bool Foam::sentinel::writeInstance(const regIOobject& io)
{
    return fileHandler().writeObject
    (
        io,
        IOstreamOption(io.time().writeFormat(), io.time().writeCompression()),
        true
    );
}


void Foam::sentinel::writeLastValid
(
    const volVectorField& U,
    const volScalarField& p,
    const surfaceScalarField& phi,
    const volScalarField* k,
    const volScalarField* omega,
    const volScalarField* nut,
    const labelList& offending,
    const label validIter
) const
{
    const word inst(lastValidName(validIter));

    auto writeCopy = [&inst](const auto& fld)
    {
        typedef std::remove_cv_t<std::remove_reference_t<decltype(fld)>> T;
        T copy
        (
            IOobject
            (
                fld.name(),
                inst,
                fld.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fld
        );
        writeInstance(copy);
    };

    writeCopy(U);
    writeCopy(p);
    writeCopy(phi);
    if (k) writeCopy(*k);
    if (omega) writeCopy(*omega);
    if (nut) writeCopy(*nut);

    volScalarField flag
    (
        IOobject
        (
            "sentinelFlag",
            inst,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );
    for (const label celli : offending)
    {
        flag[celli] = 1;
    }
    writeInstance(flag);

    // Offending cell centres, at most maxReport over all ranks
    List<pointField> allC(UPstream::nProcs());
    {
        const label n = min(offending.size(), maxReport_);
        pointField c(n);
        for (label i = 0; i < n; ++i)
        {
            c[i] = mesh_.C()[offending[i]];
        }
        allC[UPstream::myProcNo()] = c;
    }
    Pstream::gatherList(allC);

    if (UPstream::master())
    {
        label nReported = 0;
        Info<< "sentinel: offending cell centres (max " << maxReport_ << "):"
            << nl;
        forAll(allC, proci)
        {
            for (const point& pt : allC[proci])
            {
                if (nReported >= maxReport_) break;
                Info<< "    proc " << proci << "  " << pt << nl;
                ++nReported;
            }
        }
        Info<< endl;
    }
}


void Foam::sentinel::writeState(dictionary& dict) const
{
    dict.set("nRollbacks", nRollbacks_);
}


void Foam::sentinel::readState(const dictionary& dict)
{
    nRollbacks_ = dict.getOrDefault<label>("nRollbacks", 0);
    consecutive_ = 0;
}


void Foam::sentinel::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("maxRollbacks", maxRollbacks_);
    d.add("UFactor", UFactor_);
    d.add("pFactor", pFactor_);
    d.add("cflFactor", cflFactor_);
    d.add("maxReport", maxReport_);
    dict.add("sentinel", d);
}


// ************************************************************************* //
