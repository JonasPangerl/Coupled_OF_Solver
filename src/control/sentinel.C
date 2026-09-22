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

    // Turbulence fields (D-069 F3): every volScalarField the solver writes
    // except p - the model's own fields (k, omega, epsilon, nuTilda,
    // ReThetat, gammaInt, ...) and nut. sortedNames: identical order on all
    // ranks.
    if (d.found("turbulenceFields"))
    {
        turbNames_ = d.get<wordList>("turbulenceFields");
        // k, omega and nut are always included (D-056)
        for (const char* nm : {"k", "omega", "nut"})
        {
            if
            (
                mesh_.foundObject<volScalarField>(nm)
             && !turbNames_.found(word(nm))
            )
            {
                turbNames_.append(word(nm));
            }
        }
        for (const word& nm : turbNames_)
        {
            if (!mesh_.foundObject<volScalarField>(nm))
            {
                FatalIOErrorInFunction(d)
                    << "sentinel.turbulenceFields: no volScalarField "
                    << nm << " registered" << exit(FatalIOError);
            }
        }
    }
    else
    {
        DynamicList<word> names;
        for (const word& nm : mesh_.sortedNames<volScalarField>())
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(nm);
            if (nm != "p" && f.writeOpt() == IOobject::AUTO_WRITE)
            {
                names.append(nm);
            }
        }
        turbNames_.transfer(names);
    }
    Info<< "coupledFoam: sentinel turbulence fields " << flatOutput(turbNames_)
        << endl;
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
    // U and p boundary values, restored verbatim (D-069 F6)
    U0Boundary_.resize(U.boundaryField().size());
    forAll(U.boundaryField(), patchi)
    {
        U0Boundary_[patchi] = U.boundaryField()[patchi];
    }
    p0Boundary_.resize(p.boundaryField().size());
    forAll(p.boundaryField(), patchi)
    {
        p0Boundary_[patchi] = p.boundaryField()[patchi];
    }
    // All turbulence fields (k, omega, nut among them), D-069 F3
    turb0_.resize(turbNames_.size());
    turb0Boundary_.resize(turbNames_.size());
    forAll(turbNames_, i)
    {
        storeTurbulence
        (
            mesh_.lookupObject<volScalarField>(turbNames_[i]),
            turb0_[i],
            turb0Boundary_[i]
        );
    }
    stored_ = true;
}


void Foam::sentinel::storeTurbulence
(
    const volScalarField& fld,
    scalarField& internal,
    List<scalarField>& boundary
)
{
    internal = fld.primitiveField();
    boundary.resize(fld.boundaryField().size());
    forAll(fld.boundaryField(), patchi)
    {
        boundary[patchi] = fld.boundaryField()[patchi];
    }
}


void Foam::sentinel::restoreTurbulence
(
    volScalarField& fld,
    const scalarField& internal,
    const List<scalarField>& boundary
)
{
    fld.primitiveFieldRef() = internal;
    // Forced assignment (==) of the stored values: no updateCoeffs(), no
    // evaluate(). correctBoundaryConditions() would call the updateCoeffs()
    // of omegaWallFunction/epsilonWallFunction, which look up <model>:G -
    // registered only inside turbulence->correct() (D-056)
    volScalarField::Boundary& bf = fld.boundaryFieldRef();
    forAll(bf, patchi)
    {
        bf[patchi] == boundary[patchi];
    }
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

    // Turbulence fields other than k and omega (checked below with min/max)
    DynamicList<const volScalarField*> others;
    for (const word& nm : turbNames_)
    {
        const volScalarField* f = mesh_.cfindObject<volScalarField>(nm);
        if (f && f != k && f != omega)
        {
            others.append(f);
        }
    }

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

        // The other turbulence fields (epsilon, nuTilda, nut, ...): finite
        // values only (D-069 F3)
        for (const volScalarField* f : others)
        {
            if (!std::isfinite((*f)[celli]))
            {
                ++nNonFinite;
                cellBad = true;
            }
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
    // U and p boundary values verbatim as well (D-069 F6):
    // correctBoundaryConditions() would re-evaluate freestreamVelocity from
    // its own (rejected-step) patch values and inletOutlet with the restored
    // phi, i.e. give a state that is neither iteration n-1 nor the rejected
    // one. Processor patches hold the neighbour cell values of the stored
    // state (all ranks roll back together). The mixed-type coefficients are
    // refreshed by updateCoeffs() before the next assembly (F1).
    forAll(U.boundaryField(), patchi)
    {
        U.boundaryFieldRef()[patchi] == U0Boundary_[patchi];
    }
    forAll(p.boundaryField(), patchi)
    {
        p.boundaryFieldRef()[patchi] == p0Boundary_[patchi];
    }

    // All turbulence fields: verbatim, without boundary evaluation (D-056,
    // D-069 F3)
    forAll(turb0_, i)
    {
        restoreTurbulence
        (
            mesh_.lookupObjectRef<volScalarField>(turbNames_[i]),
            turb0_[i],
            turb0Boundary_[i]
        );
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
    // All turbulence fields (k, omega, nut among them), D-069 F3
    for (const word& nm : turbNames_)
    {
        writeCopy(mesh_.lookupObject<volScalarField>(nm));
    }

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
    d.add("turbulenceFields", turbNames_);
    dict.add("sentinel", d);
}


// ************************************************************************* //
