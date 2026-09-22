/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled pressure-velocity solver for OpenFOAM v2606
-------------------------------------------------------------------------------
Application
    coupledFieldCompare

Description
    Mean-field delta comparison of two runs on the same mesh (D-042
    addendum): coupledFoam vs the simpleFoam reference on the wake cases.

    Reads UMean and pMean (and UPrime2Mean where both cases have it) of this
    case at -time and of the reference case at -referenceTime (default: its
    latest time) and writes into this case at -time:
      - UMeanDelta        = UMean - UMean_ref
      - pMeanDelta        = pMean - pMean_ref
      - magUMeanDeltaRel  = |UMeanDelta| / U_inf
      - CpMeanDelta       = pMeanDelta / p_ref   (surface delta Cp on walls)
      - UPrime2MeanDelta  (only if both cases have UPrime2Mean)
    Boundary values are included (non-coupled patches: own minus reference
    patch values; coupled patches: evaluated from the delta field).

    Metrics (all reductions global) go to <case>/fieldCompare.json (master):
      - volume-weighted RMS and max of |dUMean|/U_inf and |dpMean|/p_ref
      - fraction of cells (and of the volume) with |dUMean| > 0.05 U_inf
      - per wall patch and over all walls: area-weighted RMS, mean and
        max |.| of dpMean/p_ref (the body delta Cp)
      - with UPrime2Mean: volume RMS of 0.5 tr(dUPrime2Mean)/U_inf^2

    Parallel: run with -parallel on both cases decomposed identically (same
    mesh, same decomposition); every rank reads processorN of the reference.
    The number of cells per rank must agree and the cell centres must agree
    within 1e-9 relative to max(|C|, bounding-box span); otherwise fatal.

Usage
    coupledFieldCompare -reference <caseDir> -Uinf <U> -pref <p>
        [-time <t>] [-referenceTime <t>] [-parallel]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "OFstream.H"
#include "HashTable.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Relative threshold on |dU| for the "large deviation" cell fraction
constexpr scalar largeDeviation = 0.05;

// Cell-centre agreement (relative to max(|C|, bounding-box span))
constexpr scalar centreTol = 1e-9;


//- The time directory of `runTime` that matches `t` (latest if t < 0)
instant selectTime(const Time& runTime, const scalar t, const word& what)
{
    instantList all = runTime.times();
    DynamicList<instant> times(all.size());
    for (const instant& inst : all)
    {
        scalar v = 0;
        if (readScalar(inst.name(), v))
        {
            times.append(inst);
        }
    }
    if (times.empty())
    {
        FatalErrorInFunction
            << "No time directories in " << what << " case "
            << runTime.path() << exit(FatalError);
    }
    if (t < 0)
    {
        return times.last();
    }
    for (const instant& inst : times)
    {
        if (mag(inst.value() - t) <= 1e-6*max(scalar(1), mag(t)))
        {
            return inst;
        }
    }
    FatalErrorInFunction
        << "Time " << t << " not found in " << what << " case "
        << runTime.path() << " (times: " << times << ")"
        << exit(FatalError);
    return times.last();
}


//- True if the field file exists (header check) on every rank
template<class FieldType>
bool fieldExists(const word& name, const fvMesh& mesh)
{
    IOobject io
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    return returnReduceAnd(io.typeHeaderOk<FieldType>(true));
}


template<class FieldType>
tmp<FieldType> readField(const word& name, const fvMesh& mesh)
{
    if (!fieldExists<FieldType>(name, mesh))
    {
        FatalErrorInFunction
            << "No field " << name << " at time " << mesh.time().timeName()
            << " in " << mesh.time().path() << exit(FatalError);
    }
    return tmp<FieldType>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
}


//- own - ref on `mesh` (own mesh): internal field, non-coupled patch values
//- of both fields, coupled patches evaluated
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> delta
(
    const word& name,
    const GeometricField<Type, fvPatchField, volMesh>& own,
    const GeometricField<Type, fvPatchField, volMesh>& ref
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    const fvMesh& mesh = own.mesh();

    auto td = tmp<FieldType>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>(own.dimensions(), Zero),
        fvPatchFieldBase::calculatedType()
    );
    FieldType& d = td.ref();
    d.primitiveFieldRef() = own.primitiveField() - ref.primitiveField();

    auto& bf = d.boundaryFieldRef();
    forAll(bf, patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (p.coupled() || bf[patchi].empty())
        {
            continue;
        }
        bf[patchi] ==
            own.boundaryField()[patchi] - ref.boundaryField()[patchi];
    }
    d.correctBoundaryConditions();
    return td;
}


//- Volume-weighted RMS of a cell field
scalar volRms(const scalarField& f, const scalarField& V, const scalar sumV)
{
    return Foam::sqrt(gSum(V*sqr(f))/sumV);
}


void writeNum(Ostream& os, const scalar v)
{
    if (std::isfinite(v))
    {
        os << v;
    }
    else
    {
        os << "null";
    }
}


struct wallStats
{
    word name;
    label nFaces;
    scalar area;
    scalar rms;
    scalar mean;
    scalar maxAbs;
};


wallStats wallDeltaCp
(
    const word& name,
    const scalarField& dCp,
    const scalarField& magSf
)
{
    wallStats s;
    s.name = name;
    s.nFaces = returnReduce(dCp.size(), sumOp<label>());
    s.area = gSum(magSf);
    if (s.area > 0)
    {
        s.rms = Foam::sqrt(gSum(magSf*sqr(dCp))/s.area);
        s.mean = gSum(magSf*dCp)/s.area;
        s.maxAbs = returnReduce
        (
            dCp.empty() ? scalar(0) : max(mag(dCp)),
            maxOp<scalar>()
        );
    }
    else
    {
        s.rms = s.mean = s.maxAbs = std::numeric_limits<scalar>::quiet_NaN();
    }
    return s;
}


void writeWall(Ostream& os, const wallStats& s)
{
    os  << "{\"patch\": \"" << s.name.c_str() << "\", \"nFaces\": "
        << s.nFaces << ", \"area\": ";
    writeNum(os, s.area);
    os  << ", \"rmsCpDelta\": ";
    writeNum(os, s.rms);
    os  << ", \"meanCpDelta\": ";
    writeNum(os, s.mean);
    os  << ", \"maxAbsCpDelta\": ";
    writeNum(os, s.maxAbs);
    os  << "}";
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * Main  * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mean-field delta comparison with a reference case on the same mesh"
        " (UMean, pMean, UPrime2Mean); writes delta fields and"
        " fieldCompare.json (D-042 addendum)"
    );
    argList::noFunctionObjects();
    argList::addOption
    (
        "reference",
        "dir",
        "Reference case directory (same mesh; decomposed identically when"
        " run with -parallel)"
    );
    argList::addOption("time", "t", "Time of this case (default: latest)");
    argList::addOption
    (
        "referenceTime",
        "t",
        "Time of the reference case (default: latest)"
    );
    argList::addOption("Uinf", "U", "Free-stream speed for |dU|/U_inf");
    argList::addOption("pref", "p", "Pressure scale for dp/p_ref");

    #include "setRootCase.H"

    if (!args.found("reference") || !args.found("Uinf") || !args.found("pref"))
    {
        FatalErrorInFunction
            << "-reference, -Uinf and -pref are required" << exit(FatalError);
    }
    const scalar Uinf = args.get<scalar>("Uinf");
    const scalar pref = args.get<scalar>("pref");
    if (!(Uinf > 0) || !(pref > 0))
    {
        FatalErrorInFunction
            << "-Uinf and -pref must be positive" << exit(FatalError);
    }

    // Reference database: same options, -case replaced; in parallel the
    // copied argList appends processorN like the own case
    fileName refDir = args.get<fileName>("reference");
    refDir.expand();
    if (!refDir.isAbsolute())
    {
        refDir = cwd()/refDir;
    }
    refDir.clean();
    HashTable<string> refOptions(args.options());
    refOptions.set("case", refDir);
    argList argsRef(args, refOptions, false, false, false);

    #include "createTime.H"
    Time runTimeRef(Time::controlDictName, argsRef, false, false);

    const instant tOwn =
        selectTime(runTime, args.getOrDefault<scalar>("time", -1), "own");
    const instant tRef = selectTime
    (
        runTimeRef,
        args.getOrDefault<scalar>("referenceTime", -1),
        "reference"
    );
    runTime.setTime(tOwn, 0);
    runTimeRef.setTime(tRef, 0);

    Info<< "Case      " << runTime.globalPath() << " time " << tOwn.name()
        << nl << "Reference " << runTimeRef.globalPath() << " time "
        << tRef.name() << nl << "U_inf " << Uinf << ", p_ref " << pref
        << nl << endl;

    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
    fvMesh meshRef
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTimeRef.timeName(),
            runTimeRef,
            IOobject::MUST_READ
        )
    );

    // Same mesh check: cells per rank, patches, cell centres
    {
        const bool sameCells = returnReduceAnd
        (
            mesh.nCells() == meshRef.nCells()
        );
        if (!sameCells)
        {
            FatalErrorInFunction
                << "Number of cells differs (this rank: " << mesh.nCells()
                << " vs reference " << meshRef.nCells() << "); the cases must"
                << " share the mesh and the decomposition" << exit(FatalError);
        }

        bool samePatches =
            mesh.boundaryMesh().size() == meshRef.boundaryMesh().size();
        if (samePatches)
        {
            forAll(mesh.boundaryMesh(), patchi)
            {
                const polyPatch& a = mesh.boundaryMesh()[patchi];
                const polyPatch& b = meshRef.boundaryMesh()[patchi];
                if (a.name() != b.name() || a.size() != b.size())
                {
                    samePatches = false;
                }
            }
        }
        if (!returnReduceAnd(samePatches))
        {
            FatalErrorInFunction
                << "Boundary patches differ between the cases"
                << exit(FatalError);
        }

        const vectorField& C = mesh.C().primitiveField();
        const vectorField& Cr = meshRef.C().primitiveField();
        const scalar span = mag(mesh.bounds().span());
        scalar worst = 0;
        forAll(C, celli)
        {
            const scalar scale = max(mag(C[celli]), span);
            worst = max(worst, mag(C[celli] - Cr[celli])/scale);
        }
        reduce(worst, maxOp<scalar>());
        Info<< "Max cell-centre deviation (relative) " << worst << endl;
        if (worst > centreTol)
        {
            FatalErrorInFunction
                << "Cell centres differ by " << worst << " (relative) > "
                << centreTol << ": not the same mesh/decomposition"
                << exit(FatalError);
        }
    }

    // Fields
    tmp<volVectorField> tU = readField<volVectorField>("UMean", mesh);
    tmp<volVectorField> tUr = readField<volVectorField>("UMean", meshRef);
    tmp<volScalarField> tp = readField<volScalarField>("pMean", mesh);
    tmp<volScalarField> tpr = readField<volScalarField>("pMean", meshRef);

    tmp<volVectorField> tdU = delta("UMeanDelta", tU(), tUr());
    tmp<volScalarField> tdp = delta("pMeanDelta", tp(), tpr());
    const volVectorField& dU = tdU();
    const volScalarField& dp = tdp();

    volScalarField magRel
    (
        IOobject
        (
            "magUMeanDeltaRel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(dU)/dimensionedScalar("Uinf", dU.dimensions(), Uinf)
    );
    volScalarField dCp
    (
        IOobject
        (
            "CpMeanDelta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dp/dimensionedScalar("pref", dp.dimensions(), pref)
    );

    // Volume metrics
    const scalarField& V = mesh.V().field();
    const scalar sumV = gSum(V);
    const scalarField& uRel = magRel.primitiveField();
    const scalarField pRel(mag(dCp.primitiveField()));

    const scalar rmsU = volRms(uRel, V, sumV);
    const scalar maxU = gMax(uRel);
    const scalar rmsP = volRms(pRel, V, sumV);
    const scalar maxP = gMax(pRel);

    label nLarge = 0;
    scalar vLarge = 0;
    forAll(uRel, celli)
    {
        if (uRel[celli] > largeDeviation)
        {
            ++nLarge;
            vLarge += V[celli];
        }
    }
    const label nCellsTotal = returnReduce(mesh.nCells(), sumOp<label>());
    reduce(nLarge, sumOp<label>());
    reduce(vLarge, sumOp<scalar>());

    // UPrime2Mean (optional)
    bool haveR =
        fieldExists<volSymmTensorField>("UPrime2Mean", mesh)
     && fieldExists<volSymmTensorField>("UPrime2Mean", meshRef);
    scalar rmsK = std::numeric_limits<scalar>::quiet_NaN();
    if (haveR)
    {
        tmp<volSymmTensorField> tR =
            readField<volSymmTensorField>("UPrime2Mean", mesh);
        tmp<volSymmTensorField> tRr =
            readField<volSymmTensorField>("UPrime2Mean", meshRef);
        tmp<volSymmTensorField> tdR = delta("UPrime2MeanDelta", tR(), tRr());
        const scalarField dk(0.5*tr(tdR().primitiveField())/sqr(Uinf));
        rmsK = volRms(dk, V, sumV);
        tdR().write();
    }

    // Wall patches: delta Cp
    DynamicList<wallStats> walls;
    {
        DynamicList<scalar> allCp;
        DynamicList<scalar> allA;
        forAll(mesh.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            if (!isA<wallPolyPatch>(pp))
            {
                continue;
            }
            const scalarField& cpw = dCp.boundaryField()[patchi];
            const scalarField& magSf = mesh.magSf().boundaryField()[patchi];
            walls.append(wallDeltaCp(pp.name(), cpw, magSf));
            allCp.append(cpw);
            allA.append(magSf);
        }
        // All wall faces together; a rank without wall faces still takes
        // part in the global reductions
        if (returnReduceOr(walls.size() > 0))
        {
            const scalarField cpAll(std::move(allCp));
            const scalarField aAll(std::move(allA));
            walls.append(wallDeltaCp("allWalls", cpAll, aAll));
        }
    }

    // Write the delta fields at the own time
    dU.write();
    dp.write();
    magRel.write();
    dCp.write();

    Info<< "volume RMS |dUMean|/U_inf " << rmsU << ", max " << maxU << nl
        << "volume RMS |dpMean|/p_ref " << rmsP << ", max " << maxP << nl
        << "cells with |dUMean| > " << largeDeviation << " U_inf: " << nLarge
        << " of " << nCellsTotal << nl;
    for (const wallStats& s : walls)
    {
        Info<< "wall " << s.name << ": RMS dCp " << s.rms << ", max |dCp| "
            << s.maxAbs << nl;
    }

    if (Pstream::master())
    {
        const fileName jsonFile = runTime.globalPath()/"fieldCompare.json";
        OFstream os(jsonFile);
        os.precision(15);
        os  << "{" << nl
            << "  \"case\": \"" << runTime.globalPath().c_str() << "\"," << nl
            << "  \"time\": \"" << tOwn.name().c_str() << "\"," << nl
            << "  \"reference\": \"" << runTimeRef.globalPath().c_str()
            << "\"," << nl
            << "  \"referenceTime\": \"" << tRef.name().c_str() << "\"," << nl
            << "  \"nProcs\": " << Pstream::nProcs() << "," << nl
            << "  \"Uinf\": " << Uinf << "," << nl
            << "  \"pref\": " << pref << "," << nl
            << "  \"nCells\": " << nCellsTotal << "," << nl
            << "  \"volume\": " << sumV << "," << nl
            << "  \"volRmsMagUDeltaRel\": ";
        writeNum(os, rmsU);
        os  << "," << nl << "  \"volMaxMagUDeltaRel\": ";
        writeNum(os, maxU);
        os  << "," << nl << "  \"volRmsPDeltaRel\": ";
        writeNum(os, rmsP);
        os  << "," << nl << "  \"volMaxPDeltaRel\": ";
        writeNum(os, maxP);
        os  << "," << nl << "  \"largeDeviationThreshold\": " << largeDeviation
            << "," << nl << "  \"cellFractionMagUDeltaAbove\": ";
        writeNum(os, scalar(nLarge)/max(scalar(nCellsTotal), scalar(1)));
        os  << "," << nl << "  \"volFractionMagUDeltaAbove\": ";
        writeNum(os, vLarge/sumV);
        os  << "," << nl << "  \"hasUPrime2Mean\": "
            << (haveR ? "true" : "false") << "," << nl
            << "  \"volRmsKDeltaRel\": ";
        writeNum(os, rmsK);
        os  << "," << nl << "  \"walls\": [";
        forAll(walls, i)
        {
            os  << (i ? "," : "") << nl << "    ";
            writeWall(os, walls[i]);
        }
        os  << nl << "  ]" << nl << "}" << nl;
        Info<< nl << "Wrote " << jsonFile << endl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
