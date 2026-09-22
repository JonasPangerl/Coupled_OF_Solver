/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "remediation.H"
#include "coupledDefaults.H"
#include "polyMeshTools.H"
#include "primitiveMeshTools.H"
#include "unitConversion.H"
#include "calculatedFvPatchFields.H"
#include "cellSet.H"
#include "PstreamReduceOps.H"
#include "wordRes.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::remediation::remediation
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    staticEnabled_(coupledDefaults::staticEnabled),
    nonOrthThreshold_(coupledDefaults::nonOrthThreshold),
    skewThreshold_(coupledDefaults::skewThreshold),
    volRatioThreshold_(coupledDefaults::volRatioThreshold),
    aspectThreshold_(coupledDefaults::aspectThreshold),
    staticBeta_(coupledDefaults::staticBeta),
    staticCflFactor_(coupledDefaults::staticCflFactor),
    dynamicEnabled_(coupledDefaults::dynamicEnabled),
    cU_(coupledDefaults::cU),
    cp_(coupledDefaults::cp),
    cSpike_(coupledDefaults::cSpike),
    nLayers_(coupledDefaults::nLayers),
    nQuietIters_(coupledDefaults::nQuietIters),
    dynamicCflFactor_(coupledDefaults::dynamicCflFactor),
    clipToNeighbourMean_(coupledDefaults::clipToNeighbourMean),
    warnFraction_(coupledDefaults::warnFraction),
    warnInterval_(coupledDefaults::warnInterval),
    isStatic_(mesh.nCells(), false),
    age_(mesh.nCells(), -1),
    pending_(mesh.nCells(), false),
    nStatic_(0),
    nDynamic_(0),
    dynamicVersion_(0),
    zonalEnabled_(coupledDefaults::zonalEnabled),
    zonalCfl_(),
    zonalBeta_(),
    zonalSettings_()
{
    const dictionary& r = coupledDict.subOrEmptyDict("remediation");

    const dictionary& s = r.subOrEmptyDict("static");
    staticEnabled_ = s.getOrDefault<bool>("enabled", staticEnabled_);
    nonOrthThreshold_ =
        s.getOrDefault<scalar>("nonOrthThreshold", nonOrthThreshold_);
    skewThreshold_ = s.getOrDefault<scalar>("skewThreshold", skewThreshold_);
    volRatioThreshold_ =
        s.getOrDefault<scalar>("volRatioThreshold", volRatioThreshold_);
    aspectThreshold_ =
        s.getOrDefault<scalar>("aspectThreshold", aspectThreshold_);
    staticBeta_ = s.getOrDefault<scalar>("beta", staticBeta_);
    staticCflFactor_ = s.getOrDefault<scalar>("cflFactor", staticCflFactor_);
    criteria_ = staticCriteria::settings(s);                     // C1
    criteriaBits_.resize(mesh_.nCells(), 0);

    const dictionary& d = r.subOrEmptyDict("dynamic");
    dynamicEnabled_ = d.getOrDefault<bool>("enabled", dynamicEnabled_);
    cU_ = d.getOrDefault<scalar>("cU", cU_);
    cp_ = d.getOrDefault<scalar>("cp", cp_);
    cSpike_ = d.getOrDefault<scalar>("cSpike", cSpike_);
    nLayers_ = d.getOrDefault<label>("nLayers", nLayers_);
    nQuietIters_ = d.getOrDefault<label>("nQuietIters", nQuietIters_);
    dynamicCflFactor_ = d.getOrDefault<scalar>("cflFactor", dynamicCflFactor_);
    clipToNeighbourMean_ =
        d.getOrDefault<bool>("clipToNeighbourMean", clipToNeighbourMean_);

    warnFraction_ = r.getOrDefault<scalar>("warnFraction", warnFraction_);
    warnInterval_ = r.getOrDefault<label>("warnInterval", warnInterval_);

    buildZonal(r.subOrEmptyDict("zonal"));
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::remediation::neighbourMean
(
    const volVectorField& U
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Ui = U.primitiveField();

    auto tsum = tmp<vectorField>::New(mesh_.nCells(), Zero);
    vectorField& sum = tsum.ref();
    labelField n(mesh_.nCells(), Zero);

    forAll(own, facei)
    {
        sum[own[facei]] += Ui[nei[facei]];
        sum[nei[facei]] += Ui[own[facei]];
        ++n[own[facei]];
        ++n[nei[facei]];
    }

    forAll(U.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        if (Up.coupled())
        {
            const vectorField Un(Up.patchNeighbourField());
            const labelUList& fc = Up.patch().faceCells();
            forAll(fc, pf)
            {
                sum[fc[pf]] += Un[pf];
                ++n[fc[pf]];
            }
        }
    }

    forAll(sum, celli)
    {
        // GUARD: a cell without face neighbours keeps its own value
        sum[celli] = (n[celli] > 0 ? sum[celli]/scalar(n[celli]) : Ui[celli]);
    }

    return tsum;
}


void Foam::remediation::growLayer(boolList& mark) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    boolList grown(mark);

    forAll(own, facei)
    {
        if (mark[own[facei]]) grown[nei[facei]] = true;
        if (mark[nei[facei]]) grown[own[facei]] = true;
    }

    // Across processor (and other coupled) faces
    volScalarField m
    (
        IOobject
        (
            "remediation::mark",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );
    forAll(mark, celli)
    {
        m[celli] = mark[celli] ? 1 : 0;
    }
    m.correctBoundaryConditions();

    forAll(m.boundaryField(), patchi)
    {
        const fvPatchScalarField& mp = m.boundaryField()[patchi];
        if (mp.coupled())
        {
            const scalarField mn(mp.patchNeighbourField());
            const labelUList& fc = mp.patch().faceCells();
            forAll(fc, pf)
            {
                if (mn[pf] > 0) grown[fc[pf]] = true;
            }
        }
    }

    mark.transfer(grown);
}


void Foam::remediation::readZonalFactors
(
    const dictionary& entryDict,
    const dictionary& zonalDict,
    scalar& cflFactor,
    scalar& beta
) const
{
    cflFactor = entryDict.getOrDefault<scalar>
    (
        "cflFactor",
        coupledDefaults::zonalCflFactor
    );
    beta = entryDict.getOrDefault<scalar>
    (
        "beta",
        coupledDefaults::zonalBeta
    );

    // Negated comparisons also reject non-finite input
    if (!(cflFactor > 0))
    {
        FatalIOErrorInFunction(zonalDict)
            << "zonal cflFactor must be > 0, got " << cflFactor
            << " in entry " << entryDict << exit(FatalIOError);
    }
    if (!(beta >= 0 && beta <= 1))
    {
        FatalIOErrorInFunction(zonalDict)
            << "zonal beta must be in [0, 1], got " << beta
            << " in entry " << entryDict << exit(FatalIOError);
    }
}


Foam::label Foam::remediation::applyZonal
(
    const boolList& mark,
    const scalar cflFactor,
    const scalar beta
)
{
    label nAffected = 0;
    forAll(mark, celli)
    {
        if (mark[celli])
        {
            zonalCfl_[celli] *= cflFactor;
            zonalBeta_[celli] *= beta;
            ++nAffected;
        }
    }
    return returnReduce(nAffected, sumOp<label>());
}


void Foam::remediation::buildZonal(const dictionary& zonalDict)
{
    zonalEnabled_ = zonalDict.getOrDefault<bool>("enabled", zonalEnabled_);

    zonalSettings_.clear();
    zonalSettings_.add("enabled", zonalEnabled_);

    if (!zonalEnabled_)
    {
        zonalCfl_.clear();
        zonalBeta_.clear();
        return;
    }

    zonalCfl_ = scalarField(mesh_.nCells(), scalar(1));
    zonalBeta_ = scalarField(mesh_.nCells(), scalar(1));

    const label nTotal = returnReduce(mesh_.nCells(), sumOp<label>());
    // GUARD: nTotal >= 1
    const scalar percentPerCell = 100.0/max(scalar(nTotal), scalar(1));

    const List<dictionary> zoneDicts
    (
        zonalDict.getOrDefault<List<dictionary>>("zones", List<dictionary>())
    );
    const List<dictionary> patchDicts
    (
        zonalDict.getOrDefault<List<dictionary>>
        (
            "patchDistance",
            List<dictionary>()
        )
    );

    Info<< "remediation: zonal factors enabled, " << zoneDicts.size()
        << " cellZone and " << patchDicts.size() << " patchDistance entries"
        << endl;

    // cellZone entries

    const cellZoneMesh& czm = mesh_.cellZones();
    List<dictionary> zoneSettings(zoneDicts.size());

    forAll(zoneDicts, entryi)
    {
        const dictionary& e = zoneDicts[entryi];
        const wordRe zoneName(e.get<wordRe>("cellZone"));

        scalar cflFactor = coupledDefaults::zonalCflFactor;
        scalar beta = coupledDefaults::zonalBeta;
        readZonalFactors(e, zonalDict, cflFactor, beta);

        // Amendment C5: built-in named sets (zonalSets.H)
        if (zonalSets::isBuiltin(zoneName))
        {
            if (zonalSets::needsStaticSet(zoneName))
            {
                zonalDeferred_.append({zoneName, cflFactor, beta});
                Info<< "remediation: zonal built-in " << zoneName
                    << ": applied once the static set is built" << endl;
            }
            else
            {
                const label nAffected = applyZonal
                (
                    zonalSets::mark(mesh_, zoneName, criteria_),
                    cflFactor,
                    beta
                );
                Info<< "remediation: zonal built-in " << zoneName << ": "
                    << nAffected << " cells ("
                    << percentPerCell*scalar(nAffected) << " %), cflFactor "
                    << cflFactor << ", beta " << beta << endl;
            }
            dictionary& s = zoneSettings[entryi];
            s.add("cellZone", zoneName);
            s.add("cflFactor", cflFactor);
            s.add("beta", beta);
            continue;
        }

        const labelList zoneIDs(czm.indices(zoneName, true));

        if (!returnReduce(!zoneIDs.empty(), orOp<bool>()))
        {
            FatalIOErrorInFunction(zonalDict)
                << "Unknown cellZone " << zoneName
                << " in zonal.zones entry " << entryi << nl
                << "Valid cellZones: " << flatOutput(czm.names()) << nl
                << "Valid cellZone groups: "
                << flatOutput(czm.groupZoneIDs().sortedToc())
                << exit(FatalIOError);
        }

        boolList mark(mesh_.nCells(), false);
        for (const label zonei : zoneIDs)
        {
            for (const label celli : czm[zonei])
            {
                mark[celli] = true;
            }
        }

        const label nAffected = applyZonal(mark, cflFactor, beta);

        Info<< "remediation: zonal cellZone " << zoneName << ": "
            << nAffected << " cells (" << percentPerCell*scalar(nAffected)
            << " %), cflFactor " << cflFactor << ", beta " << beta << endl;

        dictionary& s = zoneSettings[entryi];
        s.add("cellZone", zoneName);
        s.add("cflFactor", cflFactor);
        s.add("beta", beta);
    }

    // patchDistance entries

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Processor patches come last; they are decomposition artefacts and
    // never seed a layer (layers still grow across them)
    const label nNonProc = pbm.nNonProcessor();
    const wordList allPatchNames(pbm.names());
    const wordList validPatchNames(SubList<word>(allPatchNames, nNonProc));

    List<dictionary> patchSettings(patchDicts.size());

    forAll(patchDicts, entryi)
    {
        const dictionary& e = patchDicts[entryi];
        const wordRes patchNames(e.get<wordRes>("patches"));
        const label nLayers = e.getOrDefault<label>
        (
            "nLayers",
            coupledDefaults::zonalNLayers
        );

        if (nLayers < 1)
        {
            FatalIOErrorInFunction(zonalDict)
                << "zonal patchDistance nLayers must be >= 1, got " << nLayers
                << " in entry " << entryi << exit(FatalIOError);
        }

        scalar cflFactor = coupledDefaults::zonalCflFactor;
        scalar beta = coupledDefaults::zonalBeta;
        readZonalFactors(e, zonalDict, cflFactor, beta);

        boolList selected(pbm.size(), false);

        for (const wordRe& patchName : patchNames)
        {
            bool found = false;
            for (const label patchi : pbm.indices(patchName, true))
            {
                if (patchi < nNonProc)
                {
                    selected[patchi] = true;
                    found = true;
                }
            }

            if (!returnReduce(found, orOp<bool>()))
            {
                FatalIOErrorInFunction(zonalDict)
                    << "Unknown patch " << patchName
                    << " in zonal.patchDistance entry " << entryi << nl
                    << "Valid patches: " << flatOutput(validPatchNames) << nl
                    << "Valid patch groups: "
                    << flatOutput(pbm.groupPatchIDs().sortedToc())
                    << exit(FatalIOError);
            }
        }

        // Layer 1: cells adjacent to the patch faces
        boolList mark(mesh_.nCells(), false);
        forAll(selected, patchi)
        {
            if (selected[patchi])
            {
                for (const label celli : pbm[patchi].faceCells())
                {
                    mark[celli] = true;
                }
            }
        }

        // Layer k+1: face neighbours of layer k (collective)
        for (label layer = 1; layer < nLayers; ++layer)
        {
            growLayer(mark);
        }

        const label nAffected = applyZonal(mark, cflFactor, beta);

        Info<< "remediation: zonal patches " << flatOutput(patchNames)
            << " nLayers " << nLayers << ": " << nAffected << " cells ("
            << percentPerCell*scalar(nAffected) << " %), cflFactor "
            << cflFactor << ", beta " << beta << endl;

        dictionary& s = patchSettings[entryi];
        s.add("patches", patchNames);
        s.add("nLayers", nLayers);
        s.add("cflFactor", cflFactor);
        s.add("beta", beta);
    }

    zonalSettings_.add("zones", zoneSettings);
    zonalSettings_.add("patchDistance", patchSettings);
}


void Foam::remediation::applyDeferredZonal()
{
    if (!zonalEnabled_)
    {
        return;
    }

    const label nTotal = returnReduce(mesh_.nCells(), sumOp<label>());
    // GUARD: nTotal >= 1
    const scalar percentPerCell = 100.0/max(scalar(nTotal), scalar(1));

    for (const zonalSets::deferredEntry& e : zonalDeferred_)
    {
        const label nAffected = applyZonal(isStatic_, e.cflFactor, e.beta);
        Info<< "remediation: zonal built-in " << e.name << ": " << nAffected
            << " cells (" << percentPerCell*scalar(nAffected)
            << " %), cflFactor " << e.cflFactor << ", beta " << e.beta
            << endl;
    }

    // Applied once (the mesh is static)
    zonalDeferred_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::remediation::buildStatic()
{
    isStatic_ = false;
    nStatic_ = 0;
    criteriaBits_ = 0;

    if (!staticEnabled_)
    {
        Info<< "remediation: static set disabled" << endl;
        applyDeferredZonal();                                    // C5
        return;
    }

    const tmp<scalarField> tortho = polyMeshTools::faceOrthogonality
    (
        mesh_,
        mesh_.faceAreas(),
        mesh_.cellCentres()
    );
    const tmp<scalarField> tskew = polyMeshTools::faceSkewness
    (
        mesh_,
        mesh_.points(),
        mesh_.faceCentres(),
        mesh_.faceAreas(),
        mesh_.cellCentres()
    );
    const tmp<scalarField> tvolR = polyMeshTools::volRatio
    (
        mesh_,
        mesh_.cellVolumes()
    );

    scalarField openness(mesh_.nCells(), Zero);
    scalarField aratio(mesh_.nCells(), Zero);
    primitiveMeshTools::cellClosedness
    (
        mesh_,
        mesh_.geometricD(),
        mesh_.faceAreas(),
        mesh_.cellVolumes(),
        openness,
        aratio
    );

    const scalarField& ortho = tortho();
    const scalarField& skew = tskew();
    const scalarField& volR = tvolR();

    const scalar cosThreshold = std::cos(degToRad(nonOrthThreshold_));
    // volRatio is min/max <= 1; GUARD: threshold > 0
    const scalar volRatioMin = 1.0/max(volRatioThreshold_, VSMALL);

    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();

    label nOrtho = 0, nSkew = 0, nVol = 0, nAspect = 0;

    forAll(ortho, facei)
    {
        const bool bOrtho = ortho[facei] < cosThreshold;
        const bool bSkew = skew[facei] > skewThreshold_;
        const bool bVol = volR[facei] < volRatioMin;

        if (bOrtho || bSkew || bVol)
        {
            isStatic_[own[facei]] = true;
            if (facei < mesh_.nInternalFaces())
            {
                isStatic_[nei[facei]] = true;
            }
        }
        nOrtho += bOrtho;
        nSkew += bSkew;
        nVol += bVol;
    }

    forAll(aratio, celli)
    {
        if (aratio[celli] > aspectThreshold_)
        {
            isStatic_[celli] = true;
            ++nAspect;
        }
    }

    // Amendment C1: topological criteria (staticCriteria.H), OR-ed
    const staticCriteria::counts nC1 =
        staticCriteria::apply(mesh_, criteria_, criteriaBits_, isStatic_);

    forAll(isStatic_, celli)
    {
        nStatic_ += isStatic_[celli];
    }

    reduce(nOrtho, sumOp<label>());
    reduce(nSkew, sumOp<label>());
    reduce(nVol, sumOp<label>());
    reduce(nAspect, sumOp<label>());
    reduce(nStatic_, sumOp<label>());

    Info<< "remediation: static set " << nStatic_ << " cells ("
        << 100.0*scalar(nStatic_)
          /max(scalar(returnReduce(mesh_.nCells(), sumOp<label>())), scalar(1))
        << " %); faces nonOrth>" << nonOrthThreshold_ << "deg: " << nOrtho
        << ", skew>" << skewThreshold_ << ": " << nSkew
        << ", volRatio>" << volRatioThreshold_ << ": " << nVol
        << "; cells aspect>" << aspectThreshold_ << ": " << nAspect << endl;
    staticCriteria::report
    (
        criteria_,
        nC1,
        returnReduce(mesh_.nCells(), sumOp<label>())
    );

    applyDeferredZonal();                                        // C5
}


void Foam::remediation::markDynamic(const labelUList& cells)
{
    for (const label celli : cells)
    {
        pending_[celli] = true;
    }
}


void Foam::remediation::updateDynamic
(
    const volVectorField& U,
    const volScalarField& p,
    const volScalarField* kPtr,
    const volScalarField* omegaPtr,
    const scalar Uref,
    const scalar pref,
    const label iter
)
{
    if (!dynamicEnabled_)
    {
        pending_ = false;
        age_ = -1;
        nDynamic_ = 0;
        return;
    }

    boolList mark(pending_);

    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();
    const tmp<vectorField> tUbar = neighbourMean(U);
    const vectorField& Ubar = tUbar();

    forAll(mark, celli)
    {
        bool m = mark[celli];

        const vector& Uc = Ui[celli];
        const bool finiteU =
            std::isfinite(Uc.x()) && std::isfinite(Uc.y())
         && std::isfinite(Uc.z());
        const bool finiteP = std::isfinite(pi[celli]);
        bool finiteT = true;
        if (kPtr) finiteT = finiteT && std::isfinite((*kPtr)[celli]);
        if (omegaPtr) finiteT = finiteT && std::isfinite((*omegaPtr)[celli]);

        if (!finiteU || !finiteP || !finiteT)
        {
            m = true;
        }
        else
        {
            m = m
             || mag(Uc) > cU_*Uref
             || pi[celli] < -cp_*pref
             || pi[celli] > cp_*pref
             || mag(Uc - Ubar[celli]) > cSpike_*Uref;
        }
        mark[celli] = m;
    }

    for (label layer = 0; layer < nLayers_; ++layer)
    {
        growLayer(mark);
    }

    // Hysteresis
    nDynamic_ = 0;
    bool changed = false;
    forAll(age_, celli)
    {
        const bool wasIn = (age_[celli] >= 0);
        if (mark[celli])
        {
            age_[celli] = 0;
        }
        else if (age_[celli] >= 0)
        {
            ++age_[celli];
            if (age_[celli] >= nQuietIters_)
            {
                age_[celli] = -1;
            }
        }
        const bool isIn = (age_[celli] >= 0);
        changed = changed || (isIn != wasIn);
        nDynamic_ += isIn;
    }

    pending_ = false;

    // Membership changed anywhere: new version (a size comparison misses
    // cells that enter and leave in the same update)
    if (returnReduceOr(changed))
    {
        ++dynamicVersion_;
    }

    label nUnion = 0;
    forAll(age_, celli)
    {
        nUnion += (age_[celli] >= 0 || isStatic_[celli]);
    }

    reduce(nDynamic_, sumOp<label>());
    reduce(nUnion, sumOp<label>());

    const label nTotal = returnReduce(mesh_.nCells(), sumOp<label>());
    if
    (
        scalar(nUnion) > warnFraction_*scalar(nTotal)
     && warnInterval_ > 0
     && iter % warnInterval_ == 0
    )
    {
        WarningInFunction
            << "Remediation sets cover " << nUnion << " of " << nTotal
            << " cells (> " << 100*warnFraction_ << " %) at iteration "
            << iter << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::remediation::beta
(
    const scalar betaGlobal
) const
{
    auto tb = tmp<scalarField>::New(mesh_.nCells(), betaGlobal);
    scalarField& b = tb.ref();

    forAll(b, celli)
    {
        if (isStatic_[celli])
        {
            b[celli] = min(b[celli], staticBeta_);
        }
        if (age_[celli] >= 0)
        {
            b[celli] = 0;
        }
    }

    if (zonalEnabled_)
    {
        // Multipliers in [0, 1]: dynamic cells stay at 0
        b *= zonalBeta_;
    }
    return tb;
}


Foam::tmp<Foam::scalarField> Foam::remediation::cflFactor() const
{
    auto tf = tmp<scalarField>::New(mesh_.nCells(), scalar(1));
    scalarField& f = tf.ref();

    forAll(f, celli)
    {
        if (isStatic_[celli])
        {
            f[celli] = min(f[celli], staticCflFactor_);
        }
        if (age_[celli] >= 0)
        {
            f[celli] = min(f[celli], dynamicCflFactor_);
        }
    }

    if (zonalEnabled_)
    {
        f *= zonalCfl_;
    }
    return tf;
}


Foam::label Foam::remediation::clipIncrement
(
    blockScalarUList& dx,
    const volVectorField& U,
    const scalar omega,
    const scalar Uref
) const
{
    if (!dynamicEnabled_ || !clipToNeighbourMean_)
    {
        return 0;
    }

    const tmp<vectorField> tUbar = neighbourMean(U);
    const vectorField& Ubar = tUbar();
    const vectorField& Ui = U.primitiveField();
    const scalar lim = cSpike_*Uref;
    // GUARD: omega > 0 (>= omegaMin)
    const scalar rOmega = 1.0/max(omega, VSMALL);

    label nClipped = 0;

    forAll(age_, celli)
    {
        if (age_[celli] < 0)
        {
            continue;
        }

        blockScalar* d = dx.data() + celli*blockDim;
        const vector dU(d[0], d[1], d[2]);
        const vector Unew = Ui[celli] + omega*dU;
        const vector e = Unew - Ubar[celli];
        const scalar me = mag(e);

        if (me > lim)
        {
            // GUARD: me > lim >= 0
            const vector Uclip = Ubar[celli] + e*(lim/me);
            const vector dUclip = (Uclip - Ui[celli])*rOmega;
            for (label c = 0; c < blockP; ++c)
            {
                d[c] = narrow(dUclip[c]);
            }
            ++nClipped;
        }
    }

    return nClipped;
}


void Foam::remediation::write() const
{
    write(mesh_.time().timeName());
}


void Foam::remediation::write(const word& instance) const
{
    volScalarField flag
    (
        IOobject
        (
            "remediationFlag",
            instance,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );

    labelHashSet stat, dyn;
    forAll(flag, celli)
    {
        const bool s = isStatic_[celli];
        const bool d = age_[celli] >= 0;
        // Bits OR-ed (C1): 1 static, 2 dynamic, 4/8/16 static criteria
        flag[celli] = scalar
        (
            (s ? label(staticCriteria::bitStatic) : 0)
          | (d ? label(staticCriteria::bitDynamic) : 0)
          | criteriaBits_[celli]
        );
        if (s) stat.insert(celli);
        if (d) dyn.insert(celli);
    }
    flag.correctBoundaryConditions();
    flag.write();

    cellSet cs(mesh_, "remediationStatic", stat);
    cs.instance() = instance;
    cs.write();

    cellSet cd(mesh_, "remediationDynamic", dyn);
    cd.instance() = instance;
    cd.write();

    if (zonalEnabled_)
    {
        volScalarField zf
        (
            IOobject
            (
                "zonalCflFactor",
                instance,
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            calculatedFvPatchScalarField::typeName
        );
        zf.primitiveFieldRef() = zonalCfl_;
        zf.correctBoundaryConditions();
        zf.write();
    }
}


void Foam::remediation::writeState(dictionary& dict) const
{
    DynamicList<label> cells, ages;
    forAll(age_, celli)
    {
        if (age_[celli] >= 0)
        {
            cells.append(celli);
            ages.append(age_[celli]);
        }
    }
    dict.set("dynamicSet", labelList(cells));
    dict.set("dynamicSetAge", labelList(ages));
}


void Foam::remediation::readState(const dictionary& dict)
{
    age_ = -1;
    const labelList cells(dict.getOrDefault<labelList>("dynamicSet", labelList()));
    const labelList ages(dict.getOrDefault<labelList>("dynamicSetAge", labelList()));

    if (cells.size() != ages.size())
    {
        FatalIOErrorInFunction(dict)
            << "dynamicSet and dynamicSetAge differ in size"
            << exit(FatalIOError);
    }

    forAll(cells, i)
    {
        if (cells[i] < 0 || cells[i] >= mesh_.nCells())
        {
            FatalIOErrorInFunction(dict)
                << "dynamicSet cell " << cells[i] << " out of range (restart"
                << " with a different decomposition?)" << exit(FatalIOError);
        }
        age_[cells[i]] = ages[i];
    }

    nDynamic_ = returnReduce(cells.size(), sumOp<label>());
}


void Foam::remediation::writeSettings(dictionary& dict) const
{
    dictionary s;
    s.add("enabled", staticEnabled_);
    s.add("nonOrthThreshold", nonOrthThreshold_);
    s.add("skewThreshold", skewThreshold_);
    s.add("volRatioThreshold", volRatioThreshold_);
    s.add("aspectThreshold", aspectThreshold_);
    s.add("beta", staticBeta_);
    s.add("cflFactor", staticCflFactor_);
    criteria_.write(s);                                          // C1

    dictionary d;
    d.add("enabled", dynamicEnabled_);
    d.add("cU", cU_);
    d.add("cp", cp_);
    d.add("cSpike", cSpike_);
    d.add("nLayers", nLayers_);
    d.add("nQuietIters", nQuietIters_);
    d.add("cflFactor", dynamicCflFactor_);
    d.add("clipToNeighbourMean", clipToNeighbourMean_);

    dictionary r;
    r.add("static", s);
    r.add("dynamic", d);
    r.add("zonal", zonalSettings_);
    r.add("warnFraction", warnFraction_);
    r.add("warnInterval", warnInterval_);
    dict.add("remediation", r);
}


// ************************************************************************* //
