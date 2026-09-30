/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "remediation.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "polyMeshTools.H"
#include "primitiveMeshTools.H"
#include "unitConversion.H"
#include "calculatedFvPatchFields.H"
#include "cellSet.H"
#include "PstreamReduceOps.H"
#include "wordRes.H"
#include "sentinel.H"
#include <cctype>
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::remediation::remediation
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    cat_(),
    nonOrthThreshold_(coupledDefaults::nonOrthThreshold),
    skewThreshold_(coupledDefaults::skewThreshold),
    volRatioThreshold_(coupledDefaults::volRatioThreshold),
    aspectThreshold_(coupledDefaults::aspectThreshold),
    closednessThreshold_(scalar(coupledDefaults::closednessThreshold)),
    severeNonOrthThreshold_(coupledDefaults::severeNonOrthThreshold),
    severeSkewThreshold_(coupledDefaults::severeSkewThreshold),
    procNLayers_(coupledDefaults::processorNLayers),
    wallPatches_(),
    wallNLayers_(coupledDefaults::wallNLayers),
    criteria_(),
    staticBits_(mesh.nCells(), 0),
    staticBeta_(mesh.nCells(), scalar(1)),
    staticCfl_(mesh.nCells(), scalar(1)),
    gradLimited_(mesh.nCells(), false),
    nonOrthLimited_(mesh.nCells(), false),
    nCat_(label(0)),
    dynamicEnabled_(coupledDefaults::dynamicEnabled),
    cU_(coupledDefaults::cU),
    cp_(coupledDefaults::cp),
    cSpike_(coupledDefaults::cSpike),
    nLayers_(coupledDefaults::nLayers),
    nQuietIters_(coupledDefaults::nQuietIters),
    stickyAfter_(coupledDefaults::dynamicStickyAfter),
    releaseIters_(coupledDefaults::dynamicReleaseIters),
    dynamicCflFactor_(scalar(coupledDefaults::dynamicCflFactor)),
    clipToNeighbourMean_(coupledDefaults::clipToNeighbourMean),
    warnFraction_(scalar(coupledDefaults::warnFraction)),
    warnInterval_(coupledDefaults::warnInterval),
    isStatic_(mesh.nCells(), false),
    age_(mesh.nCells(), -1),
    entries_(mesh.nCells(), 0),
    nSticky_(0),
    rampLeft_(mesh.nCells(), 0),
    nRamping_(0),
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

    readCategories(r);                                           // D-066

    const dictionary& d = r.subOrEmptyDict("dynamic");
    dynamicEnabled_ = d.getOrDefault<bool>("enabled", dynamicEnabled_);
    cU_ = d.getOrDefault<scalar>("cU", cU_);
    cp_ = d.getOrDefault<scalar>("cp", cp_);
    cSpike_ = d.getOrDefault<scalar>("cSpike", cSpike_);
    nLayers_ = d.getOrDefault<label>("nLayers", nLayers_);
    nQuietIters_ = d.getOrDefault<label>("nQuietIters", nQuietIters_);
    stickyAfter_ = d.getOrDefault<label>("stickyAfter", stickyAfter_);
    if (stickyAfter_ < 0)
    {
        FatalIOErrorInFunction(d)
            << "dynamic.stickyAfter must be >= 0, got " << stickyAfter_
            << exit(FatalIOError);
    }
    releaseIters_ = d.getOrDefault<label>("releaseIters", releaseIters_);
    if (releaseIters_ < 0)
    {
        FatalIOErrorInFunction(d)
            << "dynamic.releaseIters must be >= 0, got " << releaseIters_
            << exit(FatalIOError);
    }
    dynamicCflFactor_ = d.getOrDefault<scalar>("cflFactor", dynamicCflFactor_);
    clipToNeighbourMean_ =
        d.getOrDefault<bool>("clipToNeighbourMean", clipToNeighbourMean_);

    warnFraction_ = r.getOrDefault<scalar>("warnFraction", warnFraction_);
    warnInterval_ = r.getOrDefault<label>("warnInterval", warnInterval_);

    buildZonal(r.subOrEmptyDict("zonal"));
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- Entry of a category dictionary; falls back to the legacy
//  remediation.static entry legacyKey (pre-D-066; empty: none), then to the
//  default. Records the legacy keys used. Both given with different values
//  is an error (an old-style override must not be ignored silently).
template<class T>
T categoryEntry
(
    const Foam::dictionary& catDict,
    const Foam::dictionary& legacyDict,
    const Foam::word& key,
    const Foam::word& legacyKey,
    const T& def,
    Foam::DynamicList<Foam::word>& legacyUsed
)
{
    const bool inLegacy = !legacyKey.empty() && legacyDict.found(legacyKey);
    if (inLegacy)
    {
        legacyUsed.push_uniq(legacyKey);
    }
    if (catDict.found(key))
    {
        const T v(catDict.get<T>(key));
        if (inLegacy && !(legacyDict.get<T>(legacyKey) == v))
        {
            FatalIOErrorInFunction(catDict)
                << "remediation." << catDict.dictName() << '.' << key
                << " and the deprecated remediation.static." << legacyKey
                << " are both given, with different values (D-066)"
                << Foam::exit(Foam::FatalIOError);
        }
        return v;
    }
    if (inLegacy)
    {
        return legacyDict.get<T>(legacyKey);
    }
    return def;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const char* Foam::remediation::categoryName(const label cat)
{
    static const char* names[nCategories] =
    {
        "meshQuality",
        "badMesh",
        "processor",
        "wall"
    };
    return names[cat];
}


Foam::label Foam::remediation::categoryBit(const label cat)
{
    static const label bits[nCategories] =
    {
        staticCriteria::bitMeshQuality,
        staticCriteria::bitBadMesh,
        staticCriteria::bitProcessor,
        staticCriteria::bitWall
    };
    return bits[cat];
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::remediation::readCategories(const dictionary& r)
{
    if (r.found("mild"))
    {
        FatalIOErrorInFunction(r)
            << "remediation.mild (pre-release D-061 tier) was replaced by the"
            << " categories of D-066: remediation.meshQuality, badMesh,"
            << " processor, wall" << exit(FatalIOError);
    }

    // Legacy remediation.static (pre-D-066)
    const dictionary& s = r.subOrEmptyDict("static");
    DynamicList<word> used;

    const dictionary& mq = r.subOrEmptyDict(categoryName(meshQuality));
    const dictionary& bm = r.subOrEmptyDict(categoryName(badMesh));
    const dictionary& pr = r.subOrEmptyDict(categoryName(processor));
    const dictionary& wa = r.subOrEmptyDict(categoryName(wall));

    // Switches and treatment. Legacy: static.enabled switches all four,
    // static.beta/cflFactor the full-treatment categories
    const FixedList<const dictionary*, nCategories> dicts({&mq, &bm, &pr, &wa});
    const FixedList<bool, nCategories> defEnabled
    ({
        coupledDefaults::meshQualityEnabled,
        coupledDefaults::badMeshEnabled,
        coupledDefaults::processorEnabled,
        coupledDefaults::wallEnabled
    });
    const FixedList<scalar, nCategories> defBeta
    ({
        coupledDefaults::meshQualityBeta,
        coupledDefaults::badMeshBeta,
        coupledDefaults::processorBeta,
        coupledDefaults::wallBeta
    });
    const FixedList<scalar, nCategories> defCfl
    ({
        coupledDefaults::meshQualityCflFactor,
        coupledDefaults::badMeshCflFactor,
        coupledDefaults::processorCflFactor,
        coupledDefaults::wallCflFactor
    });
    const FixedList<bool, nCategories> defGrad
    ({
        coupledDefaults::meshQualityGradLimiter,
        coupledDefaults::badMeshGradLimiter,
        coupledDefaults::processorGradLimiter,
        coupledDefaults::wallGradLimiter
    });
    const FixedList<bool, nCategories> defNonOrth
    ({
        coupledDefaults::meshQualityNonOrthLimiter,
        coupledDefaults::badMeshNonOrthLimiter,
        coupledDefaults::processorNonOrthLimiter,
        coupledDefaults::wallNonOrthLimiter
    });

    for (label c = 0; c < nCategories; ++c)
    {
        const dictionary& cd = *dicts[c];
        const bool full = (c == meshQuality || c == badMesh);
        categorySettings& cs = cat_[c];
        cs.enabled = categoryEntry<bool>
        (
            cd, s, "enabled", "enabled", defEnabled[c], used
        );
        cs.beta = categoryEntry<scalar>
        (
            cd, s, "beta", full ? "beta" : "", defBeta[c], used
        );
        cs.cflFactor = categoryEntry<scalar>
        (
            cd, s, "cflFactor", full ? "cflFactor" : "", defCfl[c], used
        );
        cs.gradLimiter = cd.getOrDefault<bool>("gradLimiter", defGrad[c]);
        cs.nonOrthLimiter =
            cd.getOrDefault<bool>("nonOrthLimiter", defNonOrth[c]);

        // Negated comparisons also reject non-finite input
        if (!(cs.beta >= 0 && cs.beta <= 1))
        {
            FatalIOErrorInFunction(r)
                << "remediation." << categoryName(c) << ".beta must be in"
                << " [0, 1], got " << cs.beta << exit(FatalIOError);
        }
        if (!(cs.cflFactor > 0 && cs.cflFactor <= 1))
        {
            FatalIOErrorInFunction(r)
                << "remediation." << categoryName(c) << ".cflFactor must be"
                << " in (0, 1], got " << cs.cflFactor << exit(FatalIOError);
        }
    }

    // meshQuality criteria (8.1)
    nonOrthThreshold_ = categoryEntry<scalar>
    (
        mq, s, "nonOrthThreshold", "nonOrthThreshold", nonOrthThreshold_, used
    );
    skewThreshold_ = categoryEntry<scalar>
    (
        mq, s, "skewThreshold", "skewThreshold", skewThreshold_, used
    );
    volRatioThreshold_ = categoryEntry<scalar>
    (
        mq, s, "volRatioThreshold", "volRatioThreshold", volRatioThreshold_,
        used
    );
    aspectThreshold_ = categoryEntry<scalar>
    (
        mq, s, "aspectThreshold", "aspectThreshold", aspectThreshold_, used
    );

    // badMesh criteria
    criteria_.volumeJump = categoryEntry<bool>
    (
        bm, s, "volumeJump", "volumeJump", criteria_.volumeJump, used
    );
    criteria_.volJumpThreshold = categoryEntry<scalar>
    (
        bm, s, "volJumpThreshold", "volJumpThreshold",
        criteria_.volJumpThreshold, used
    );
    closednessThreshold_ =
        bm.getOrDefault<scalar>("closednessThreshold", closednessThreshold_);
    severeNonOrthThreshold_ =
        bm.getOrDefault<scalar>("nonOrthThreshold", severeNonOrthThreshold_);
    severeSkewThreshold_ =
        bm.getOrDefault<scalar>("skewThreshold", severeSkewThreshold_);
    if
    (
        !(closednessThreshold_ >= 0)
     || !(severeNonOrthThreshold_ >= 0 && severeNonOrthThreshold_ <= 180)
     || !(severeSkewThreshold_ >= 0)
    )
    {
        FatalIOErrorInFunction(r)
            << "remediation.badMesh: closednessThreshold ("
            << closednessThreshold_ << ") and skewThreshold ("
            << severeSkewThreshold_ << ") must be >= 0, nonOrthThreshold ("
            << severeNonOrthThreshold_ << ") in [0, 180] (0 = off)"
            << exit(FatalIOError);
    }

    // processor criteria
    criteria_.procAMI = categoryEntry<bool>
    (
        pr, s, "procAMI", "procAMI", criteria_.procAMI, used
    );
    procNLayers_ = pr.getOrDefault<label>("nLayers", procNLayers_);
    if (procNLayers_ < 0)
    {
        FatalIOErrorInFunction(r)
            << "remediation.processor.nLayers must be >= 0, got "
            << procNLayers_ << exit(FatalIOError);
    }

    // wall criteria
    criteria_.wallStarved = categoryEntry<bool>
    (
        wa, s, "wallStarved", "wallStarved", criteria_.wallStarved, used
    );
    criteria_.maxWallInternalFaces = categoryEntry<label>
    (
        wa, s, "maxWallInternalFaces", "maxWallInternalFaces",
        criteria_.maxWallInternalFaces, used
    );
    wallPatches_ = wa.getOrDefault<wordRes>("patches", wordRes());
    wallNLayers_ = wa.getOrDefault<label>("nLayers", wallNLayers_);
    if (wallNLayers_ < 1)
    {
        FatalIOErrorInFunction(r)
            << "remediation.wall.nLayers must be >= 1, got " << wallNLayers_
            << exit(FatalIOError);
    }

    criteria_.check(r);

    // The legacy non-orthogonal limiter coefficient is read by the
    // assembler (staticCriteria::limitedNonOrthCoeff); listed here
    if (s.found("nonOrthLimiter") && !r.found("limitedNonOrthCoeff"))
    {
        used.push_uniq("nonOrthLimiter");
    }

    if (r.isDict("static"))
    {
        Info<< "remediation: NOTE remediation.static is deprecated (D-066);"
            << " mapped onto the categories: " << flatOutput(used);
        if (used.found("wallStarved"))
        {
            Info<< " (static.wallStarved now gets the mild wall treatment,"
                << " D-061)";
        }
        Info<< endl;
    }
}


Foam::boolList Foam::remediation::patchLayers
(
    const wordRes& patchNames,
    const label nLayers,
    const dictionary& errDict,
    const string& context
) const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Processor patches come last; they are decomposition artefacts and
    // never seed a layer (layers still grow across them)
    const label nNonProc = pbm.nNonProcessor();

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
            const wordList allPatchNames(pbm.names());
            FatalIOErrorInFunction(errDict)
                << "Unknown patch " << patchName << " in " << context << nl
                << "Valid patches: "
                << flatOutput(SubList<word>(allPatchNames, nNonProc)) << nl
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

    return mark;
}


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
    const scalar percentPerCell = scalar(100.0/max(scalar(nTotal), scalar(1)));

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

        const boolList mark
        (
            patchLayers
            (
                patchNames,
                nLayers,
                zonalDict,
                "zonal.patchDistance entry " + Foam::name(entryi)
            )
        );

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
    const scalar percentPerCell = scalar(100.0/max(scalar(nTotal), scalar(1)));

    for (const zonalSets::deferredEntry& e : zonalDeferred_)
    {
        const label bit = zonalSets::staticSetBit(e.name);
        boolList mark(mesh_.nCells(), false);
        forAll(mark, celli)
        {
            mark[celli] = (staticBits_[celli] & bit) != 0;
        }
        const label nAffected = applyZonal(mark, e.cflFactor, e.beta);
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
    nCat_ = label(0);
    staticBits_ = 0;
    staticBeta_ = scalar(1);
    staticCfl_ = scalar(1);
    gradLimited_ = false;
    nonOrthLimited_ = false;

    const label nTotal = returnReduce(mesh_.nCells(), sumOp<label>());
    // GUARD: nTotal >= 1
    const scalar percentPerCell = scalar(100.0/max(scalar(nTotal), scalar(1)));

    bool anyEnabled = false;
    for (const categorySettings& cs : cat_)
    {
        anyEnabled = anyEnabled || cs.enabled;
    }
    if (!anyEnabled)
    {
        Info<< "remediation: static set disabled (all categories off)"
            << endl;
        applyDeferredZonal();                                    // C5
        return;
    }

    const bool mqOn = cat_[meshQuality].enabled;
    const bool bmOn = cat_[badMesh].enabled;

    // Quality criteria (meshQuality and badMesh), per-face and per-cell
    label nOrtho = 0, nSkew = 0, nVol = 0, nAspect = 0;
    label nClosed = 0, nSevere = 0;

    if (mqOn || bmOn)
    {
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
        const scalar volRatioMin =
            scalar(1.0/max(volRatioThreshold_, cfVSmall<scalar>()));

        // Severe face criteria of badMesh (0 = off)
        const bool severeOrthoOn = bmOn && severeNonOrthThreshold_ > 0;
        const bool severeSkewOn = bmOn && severeSkewThreshold_ > 0;
        const scalar cosSevere =
            std::cos(degToRad(severeNonOrthThreshold_));

        const labelUList& own = mesh_.faceOwner();
        const labelUList& nei = mesh_.faceNeighbour();

        auto markFace = [&](const label facei, const label bit)
        {
            staticBits_[own[facei]] |= bit;
            if (facei < mesh_.nInternalFaces())
            {
                staticBits_[nei[facei]] |= bit;
            }
        };

        forAll(ortho, facei)
        {
            if (mqOn)
            {
                const bool bOrtho = ortho[facei] < cosThreshold;
                const bool bSkew = skew[facei] > skewThreshold_;
                const bool bVol = volR[facei] < volRatioMin;

                if (bOrtho || bSkew || bVol)
                {
                    markFace(facei, staticCriteria::bitMeshQuality);
                }
                nOrtho += bOrtho;
                nSkew += bSkew;
                nVol += bVol;
            }
            if
            (
                (severeOrthoOn && ortho[facei] < cosSevere)
             || (severeSkewOn && skew[facei] > severeSkewThreshold_)
            )
            {
                markFace(facei, staticCriteria::bitSevereQuality);
                ++nSevere;
            }
        }

        forAll(aratio, celli)
        {
            if (mqOn && aratio[celli] > aspectThreshold_)
            {
                staticBits_[celli] |= staticCriteria::bitMeshQuality;
                ++nAspect;
            }
            if
            (
                bmOn
             && closednessThreshold_ > 0
             && openness[celli] > closednessThreshold_
            )
            {
                staticBits_[celli] |= staticCriteria::bitClosedness;
                ++nClosed;
            }
        }
    }

    // Amendment C1: topological criteria (staticCriteria.H) of the enabled
    // categories. The switches are identical on all ranks, so the
    // collectives inside are called everywhere or nowhere.
    staticCriteria::settings c1(criteria_);
    c1.volumeJump = c1.volumeJump && bmOn;
    c1.procAMI = c1.procAMI && cat_[processor].enabled;
    c1.wallStarved = c1.wallStarved && cat_[wall].enabled;
    const staticCriteria::counts nC1 =
        staticCriteria::apply(mesh_, c1, staticBits_);

    // processor: layers from the processor patches (collective)
    label nProcLayer = 0;
    if (cat_[processor].enabled && procNLayers_ > 0)
    {
        boolList mark(staticCriteria::procCells(mesh_));
        for (label layer = 1; layer < procNLayers_; ++layer)
        {
            growLayer(mark);
        }
        forAll(mark, celli)
        {
            if (mark[celli])
            {
                staticBits_[celli] |= staticCriteria::bitProcLayer;
                ++nProcLayer;
            }
        }
    }

    // wall: layers from the wall.patches patches (collective)
    label nWallLayer = 0;
    if (cat_[wall].enabled && !wallPatches_.empty())
    {
        const boolList mark
        (
            patchLayers
            (
                wallPatches_,
                wallNLayers_,
                mesh_.solutionDict(),
                "coupled.remediation.wall.patches"
            )
        );
        forAll(mark, celli)
        {
            if (mark[celli])
            {
                staticBits_[celli] |= staticCriteria::bitWallLayer;
                ++nWallLayer;
            }
        }
    }

    // Criteria -> categories (only enabled categories have criteria bits)
    const FixedList<label, nCategories> criteriaOf
    ({
        label(staticCriteria::bitMeshQuality),
        label
        (
            staticCriteria::bitVolumeJump
          | staticCriteria::bitClosedness
          | staticCriteria::bitSevereQuality
        ),
        label(staticCriteria::bitProcAMI | staticCriteria::bitProcLayer),
        label(staticCriteria::bitWallStarved | staticCriteria::bitWallLayer)
    });

    // Per-cell treatment: the most restrictive value of the cell's
    // categories (min beta, min cflFactor, OR of the limiter switches)
    forAll(staticBits_, celli)
    {
        label& bits = staticBits_[celli];
        for (label c = 0; c < nCategories; ++c)
        {
            if (cat_[c].enabled && (bits & criteriaOf[c]))
            {
                const categorySettings& cs = cat_[c];
                bits |= categoryBit(c) | staticCriteria::bitStatic;
                staticBeta_[celli] = min(staticBeta_[celli], cs.beta);
                staticCfl_[celli] = min(staticCfl_[celli], cs.cflFactor);
                gradLimited_[celli] = gradLimited_[celli] || cs.gradLimiter;
                nonOrthLimited_[celli] =
                    nonOrthLimited_[celli] || cs.nonOrthLimiter;
                ++nCat_[c];
            }
        }
        isStatic_[celli] = (bits & staticCriteria::bitStatic) != 0;
        nStatic_ += isStatic_[celli];
    }

    reduce(nOrtho, sumOp<label>());
    reduce(nSkew, sumOp<label>());
    reduce(nVol, sumOp<label>());
    reduce(nAspect, sumOp<label>());
    reduce(nClosed, sumOp<label>());
    reduce(nSevere, sumOp<label>());
    reduce(nProcLayer, sumOp<label>());
    reduce(nWallLayer, sumOp<label>());
    reduce(nStatic_, sumOp<label>());
    for (label& n : nCat_)
    {
        reduce(n, sumOp<label>());
    }

    // Once-per-run block (12.2)

    auto pct = [percentPerCell](const label n)
    {
        return percentPerCell*scalar(n);
    };
    auto item = [&](const bool on, const label n)
    {
        if (on)
        {
            Info<< n;
        }
        else
        {
            Info<< "off";
        }
    };
    auto yesNo = [](const bool b) { return b ? "yes" : "no"; };
    auto head = [&](const label c)
    {
        const categorySettings& cs = cat_[c];
        Info<< "remediation:   " << categoryName(c) << ": ";
        if (!cs.enabled)
        {
            Info<< "off" << endl;
            return false;
        }
        Info<< nCat_[c] << " cells (" << pct(nCat_[c]) << " %), beta "
            << cs.beta << ", cflFactor " << cs.cflFactor << ", gradLimiter "
            << yesNo(cs.gradLimiter) << ", nonOrthLimiter "
            << yesNo(cs.nonOrthLimiter) << "; ";
        return true;
    };

    Info<< "remediation: static set " << nStatic_ << " cells ("
        << pct(nStatic_) << " %) in the categories (D-066):" << endl;
    if (head(meshQuality))
    {
        Info<< "faces nonOrth>" << nonOrthThreshold_ << "deg: " << nOrtho
            << ", skew>" << skewThreshold_ << ": " << nSkew
            << ", volRatio>" << volRatioThreshold_ << ": " << nVol
            << "; cells aspect>" << aspectThreshold_ << ": " << nAspect
            << endl;
    }
    if (head(badMesh))
    {
        Info<< "cells volumeJump>" << c1.volJumpThreshold << ": ";
        item(c1.volumeJump, nC1.nVolumeJump);
        Info<< ", closedness>" << closednessThreshold_ << ": ";
        item(closednessThreshold_ > 0, nClosed);
        Info<< "; faces nonOrth>" << severeNonOrthThreshold_ << "deg or skew>"
            << severeSkewThreshold_ << ": ";
        item
        (
            severeNonOrthThreshold_ > 0 || severeSkewThreshold_ > 0,
            nSevere
        );
        Info<< endl;
    }
    if (head(processor))
    {
        Info<< "cells procAMI: ";
        item(c1.procAMI, nC1.nProcAMI);
        Info<< ", processor-patch layers (nLayers " << procNLayers_ << "): ";
        item(procNLayers_ > 0, nProcLayer);
        Info<< endl;
    }
    if (head(wall))
    {
        Info<< "cells wallStarved (<= " << c1.maxWallInternalFaces
            << " internal faces, 3D): ";
        item(c1.wallStarved, nC1.nWallStarved);
        Info<< ", patch layers " << flatOutput(wallPatches_) << " (nLayers "
            << wallNLayers_ << "): ";
        item(!wallPatches_.empty(), nWallLayer);
        Info<< endl;
    }

    applyDeferredZonal();                                        // C5
}


void Foam::remediation::markDynamic(const labelUList& cells)
{
    for (const label celli : cells)
    {
        pending_[celli] = true;
    }
}


Foam::label Foam::remediation::activatePending(const label extraLayers)
{
    // D-076. Before this existed, a rollback only set pending_, and pending_
    // is consumed by updateDynamic - which runs after an ACCEPTED step, and
    // not at all during a start-up with startupReference exclude. In a
    // cascade of rollbacks no step is accepted, so the cells the sentinel
    // had just named were never treated, and every retry was assembled
    // exactly like the step that had failed (F1 half-car 2026-09-30: nDyn=0
    // through five consecutive rollbacks).
    if (!dynamicEnabled_)
    {
        pending_ = false;
        return 0;
    }

    boolList mark(pending_);
    pending_ = false;

    // Collective: the layer count is the same on every rank (the rollback
    // decision and the consecutive-rollback count are global)
    const label nLayers = nLayers_ + max(extraLayers, label(0));
    for (label layer = 0; layer < nLayers; ++layer)
    {
        growLayer(mark);
    }

    label nEntered = 0;
    bool changed = false;
    nDynamic_ = 0;
    nSticky_ = 0;
    nRamping_ = 0;
    forAll(age_, celli)
    {
        if (mark[celli])
        {
            if (age_[celli] < 0)
            {
                ++entries_[celli];
                ++nEntered;
                changed = true;
            }
            if (rampLeft_[celli] > 0)
            {
                // Leaves its release ramp: beta and dt change
                changed = true;
            }
            // A member restarts its quiet count; other cells do not age,
            // no step was accepted
            age_[celli] = 0;
            rampLeft_[celli] = 0;
        }
        const bool isIn = (age_[celli] >= 0);
        const bool sticky =
            (stickyAfter_ > 0 && entries_[celli] >= stickyAfter_);
        nDynamic_ += isIn;
        nSticky_ += (isIn && sticky);
        nRamping_ += (rampLeft_[celli] > 0);
    }

    reduce(nDynamic_, sumOp<label>());
    reduce(nSticky_, sumOp<label>());
    reduce(nRamping_, sumOp<label>());
    reduce(nEntered, sumOp<label>());

    if (returnReduceOr(changed))
    {
        ++dynamicVersion_;
    }

    return nEntered;
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
        rampLeft_ = 0;
        nDynamic_ = 0;
        nRamping_ = 0;
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

    // Hysteresis; sticky cells (stickyAfter_ entries) are never released
    nDynamic_ = 0;
    nSticky_ = 0;
    nRamping_ = 0;
    bool changed = false;
    forAll(age_, celli)
    {
        const bool wasIn = (age_[celli] >= 0);
        // A running release ramp advances by one iteration (its beta and
        // dt change: a new operator, as for a membership change)
        if (rampLeft_[celli] > 0)
        {
            --rampLeft_[celli];
            changed = true;
        }
        if (mark[celli])
        {
            if (!wasIn)
            {
                ++entries_[celli];
            }
            age_[celli] = 0;
            rampLeft_[celli] = 0;
        }
        const bool sticky =
            (stickyAfter_ > 0 && entries_[celli] >= stickyAfter_);
        if (!mark[celli] && age_[celli] >= 0)
        {
            ++age_[celli];
            if (age_[celli] >= nQuietIters_ && !sticky)
            {
                age_[celli] = -1;
                rampLeft_[celli] = releaseIters_;
            }
        }
        const bool isIn = (age_[celli] >= 0);
        changed = changed || (isIn != wasIn);
        nDynamic_ += isIn;
        nSticky_ += (isIn && sticky);
        nRamping_ += (rampLeft_[celli] > 0);
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
    reduce(nSticky_, sumOp<label>());
    reduce(nRamping_, sumOp<label>());
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
            // Category treatment (D-066)
            b[celli] = min(b[celli], staticBeta_[celli]);
        }
        if (age_[celli] >= 0)
        {
            b[celli] = 0;
        }
        else if (rampLeft_[celli] > 0)
        {
            // Release ramp: 1 - j/(N+1), j = iterations left
            b[celli] *=
                1 - scalar(rampLeft_[celli])/scalar(releaseIters_ + 1);
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
            // Category treatment (D-066)
            f[celli] = min(f[celli], staticCfl_[celli]);
        }
        if (age_[celli] >= 0)
        {
            f[celli] = min(f[celli], dynamicCflFactor_);
        }
        else if (rampLeft_[celli] > 0)
        {
            // Release ramp: cflFactor^(j/(N+1)) (GUARD: cflFactor > 0)
            f[celli] = min
            (
                f[celli],
                std::pow
                (
                    max(dynamicCflFactor_, cfVSmall<scalar>()),
                    scalar(rampLeft_[celli])/scalar(releaseIters_ + 1)
                )
            );
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
    const scalar rOmega = scalar(1.0/max(omega, cfVSmall<scalar>()));

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
    FixedList<labelHashSet, nCategories> cats;
    forAll(flag, celli)
    {
        const bool s = isStatic_[celli];
        const bool d = age_[celli] >= 0;
        // Bits OR-ed (staticCriteria::flagBits): 1 static, 2 dynamic,
        // criteria and categories (C1, D-066)
        flag[celli] = scalar
        (
            (d ? label(staticCriteria::bitDynamic) : 0)
          | staticBits_[celli]
        );
        if (s) stat.insert(celli);
        if (d) dyn.insert(celli);
        for (label c = 0; c < nCategories; ++c)
        {
            if (staticBits_[celli] & categoryBit(c))
            {
                cats[c].insert(celli);
            }
        }
    }
    // Written through the file handler: regIOobject::writeObject would
    // redirect a non-time instance (<n>_lastValid) to the current time
    flag.correctBoundaryConditions();
    sentinel::writeInstance(flag);

    cellSet cs(mesh_, "remediationStatic", stat);
    cs.instance() = instance;
    sentinel::writeInstance(cs);

    cellSet cd(mesh_, "remediationDynamic", dyn);
    cd.instance() = instance;
    sentinel::writeInstance(cd);

    // One cellSet per category (D-066): remediationMeshQuality, ...
    for (label c = 0; c < nCategories; ++c)
    {
        word setName(categoryName(c));
        setName[0] = char(std::toupper(setName[0]));
        cellSet cc(mesh_, "remediation" + setName, cats[c]);
        cc.instance() = instance;
        sentinel::writeInstance(cc);
    }

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
        sentinel::writeInstance(zf);
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

    // Entry counters (sparse, local labels)
    DynamicList<label> entryCells, entryCounts;
    forAll(entries_, celli)
    {
        if (entries_[celli] > 0)
        {
            entryCells.append(celli);
            entryCounts.append(entries_[celli]);
        }
    }
    dict.set("dynamicEntryCells", labelList(entryCells));
    dict.set("dynamicEntryCount", labelList(entryCounts));

    // Release ramps (sparse, local labels)
    DynamicList<label> rampCells, rampLeft;
    forAll(rampLeft_, celli)
    {
        if (rampLeft_[celli] > 0)
        {
            rampCells.append(celli);
            rampLeft.append(rampLeft_[celli]);
        }
    }
    dict.set("dynamicRampCells", labelList(rampCells));
    dict.set("dynamicRampLeft", labelList(rampLeft));

    // Marks of a rollback or omegaMin step of the written iteration, not
    // yet applied by updateDynamic (D-069 F4)
    DynamicList<label> pendingCells;
    forAll(pending_, celli)
    {
        if (pending_[celli])
        {
            pendingCells.append(celli);
        }
    }
    dict.set("dynamicPending", labelList(pendingCells));
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

    entries_ = 0;
    const labelList entryCells
    (
        dict.getOrDefault<labelList>("dynamicEntryCells", labelList())
    );
    const labelList entryCounts
    (
        dict.getOrDefault<labelList>("dynamicEntryCount", labelList())
    );
    if (entryCells.size() != entryCounts.size())
    {
        FatalIOErrorInFunction(dict)
            << "dynamicEntryCells and dynamicEntryCount differ in size"
            << exit(FatalIOError);
    }
    forAll(entryCells, i)
    {
        if (entryCells[i] < 0 || entryCells[i] >= mesh_.nCells())
        {
            FatalIOErrorInFunction(dict)
                << "dynamicEntryCells cell " << entryCells[i]
                << " out of range (restart with a different decomposition?)"
                << exit(FatalIOError);
        }
        entries_[entryCells[i]] = entryCounts[i];
    }

    label nSticky = 0;
    forAll(age_, celli)
    {
        nSticky +=
        (
            age_[celli] >= 0
         && stickyAfter_ > 0
         && entries_[celli] >= stickyAfter_
        );
    }
    nSticky_ = returnReduce(nSticky, sumOp<label>());

    rampLeft_ = 0;
    const labelList rampCells
    (
        dict.getOrDefault<labelList>("dynamicRampCells", labelList())
    );
    const labelList rampLeft
    (
        dict.getOrDefault<labelList>("dynamicRampLeft", labelList())
    );
    if (rampCells.size() != rampLeft.size())
    {
        FatalIOErrorInFunction(dict)
            << "dynamicRampCells and dynamicRampLeft differ in size"
            << exit(FatalIOError);
    }
    forAll(rampCells, i)
    {
        if (rampCells[i] < 0 || rampCells[i] >= mesh_.nCells())
        {
            FatalIOErrorInFunction(dict)
                << "dynamicRampCells cell " << rampCells[i]
                << " out of range (restart with a different decomposition?)"
                << exit(FatalIOError);
        }
        // A ramp longer than the current releaseIters is shortened
        rampLeft_[rampCells[i]] = min(rampLeft[i], releaseIters_);
    }
    label nRamp = 0;
    forAll(rampLeft_, celli)
    {
        nRamp += (rampLeft_[celli] > 0);
    }
    nRamping_ = returnReduce(nRamp, sumOp<label>());

    // Pending marks (D-069 F4)
    pending_ = false;
    for
    (
        const label celli
      : dict.getOrDefault<labelList>("dynamicPending", labelList())
    )
    {
        if (celli < 0 || celli >= mesh_.nCells())
        {
            FatalIOErrorInFunction(dict)
                << "dynamicPending cell " << celli << " out of range"
                << " (restart with a different decomposition?)"
                << exit(FatalIOError);
        }
        pending_[celli] = true;
    }
}


void Foam::remediation::writeSettings(dictionary& dict) const
{
    // Categories (D-066)
    FixedList<dictionary, nCategories> cs;
    for (label c = 0; c < nCategories; ++c)
    {
        cs[c].add("enabled", cat_[c].enabled);
    }
    cs[meshQuality].add("nonOrthThreshold", nonOrthThreshold_);
    cs[meshQuality].add("skewThreshold", skewThreshold_);
    cs[meshQuality].add("volRatioThreshold", volRatioThreshold_);
    cs[meshQuality].add("aspectThreshold", aspectThreshold_);
    cs[badMesh].add("volumeJump", criteria_.volumeJump);
    cs[badMesh].add("volJumpThreshold", criteria_.volJumpThreshold);
    cs[badMesh].add("closednessThreshold", closednessThreshold_);
    cs[badMesh].add("nonOrthThreshold", severeNonOrthThreshold_);
    cs[badMesh].add("skewThreshold", severeSkewThreshold_);
    cs[processor].add("procAMI", criteria_.procAMI);
    cs[processor].add("nLayers", procNLayers_);
    cs[wall].add("wallStarved", criteria_.wallStarved);
    cs[wall].add("maxWallInternalFaces", criteria_.maxWallInternalFaces);
    cs[wall].add("patches", wallPatches_);
    cs[wall].add("nLayers", wallNLayers_);

    dictionary r;
    for (label c = 0; c < nCategories; ++c)
    {
        cs[c].add("beta", cat_[c].beta);
        cs[c].add("cflFactor", cat_[c].cflFactor);
        cs[c].add("gradLimiter", cat_[c].gradLimiter);
        cs[c].add("nonOrthLimiter", cat_[c].nonOrthLimiter);
        r.add(categoryName(c), cs[c]);
    }

    dictionary d;
    d.add("enabled", dynamicEnabled_);
    d.add("cU", cU_);
    d.add("cp", cp_);
    d.add("cSpike", cSpike_);
    d.add("nLayers", nLayers_);
    d.add("nQuietIters", nQuietIters_);
    d.add("stickyAfter", stickyAfter_);
    d.add("releaseIters", releaseIters_);
    d.add("cflFactor", dynamicCflFactor_);
    d.add("clipToNeighbourMean", clipToNeighbourMean_);

    r.add("dynamic", d);
    r.add("zonal", zonalSettings_);
    r.add("warnFraction", warnFraction_);
    r.add("warnInterval", warnInterval_);
    dict.add("remediation", r);
}


// ************************************************************************* //
