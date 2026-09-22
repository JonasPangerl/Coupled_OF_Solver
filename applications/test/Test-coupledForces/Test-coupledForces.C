/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-coupledForces

Description
    Unit test of coupledForces (amendment D6): on a solved state of a case
    (T3 airFoil2D in the test suite) the double-accumulated coefficients
    must equal those of the native forceCoeffs function object of the
    case's controlDict to a relative tolerance (default 1e-6) for Cd, Cl
    and Cm (= native CmPitch); the pressure/viscous parts are reported.

    Both evaluations use the same fields in the same process: p, U and the
    turbulence model read from the selected time (nut as written, no
    validate()), devReff from the model. Runs serial or -parallel (the
    double reduction of coupledForces against the native scalar one).

Usage
    Test-coupledForces [-latestTime | -time <t>] [-native <entry>]
                       [-tol <rel>] [-json <file>] [-parallel]

    -native   entry of controlDict functions holding the native
              forceCoeffs settings (default forceCoeffs)
    Settings of coupledForces: system/coupledForcesDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "forceCoeffs.H"
#include "timeSelector.H"
#include "coupledForces.H"
#include "jsonWriter.H"
#include <cmath>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "coupledForces (D6) vs native forceCoeffs on a solved state"
    );
    timeSelector::addOptions(false, false);
    argList::noFunctionObjects();
    argList::addOption("native", "entry", "native forceCoeffs entry name");
    argList::addOption("tol", "rel", "relative tolerance (default 1e-6)");
    argList::addOption("json", "file", "Write results as JSON");

    #include "setRootCase.H"
    #include "createTime.H"

    const instantList times = timeSelector::select0(runTime, args);
    if (times.empty())
    {
        FatalErrorInFunction << "no time selected" << exit(FatalError);
    }
    runTime.setTime(times.last(), times.size() - 1);

    #include "createMesh.H"

    const word nativeName(args.getOrDefault<word>("native", "forceCoeffs"));
    const doubleScalar tol(args.getOrDefault<doubleScalar>("tol", 1e-6));

    Info<< "Time = " << runTime.timeName() << nl << endl;

    volScalarField p
    (
        IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ),
        mesh
    );
    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    // --- native forceCoeffs with the case's settings
    const dictionary& fdict =
        runTime.controlDict().subDict("functions").subDict(nativeName);
    functionObjects::forceCoeffs native(nativeName, runTime, fdict);
    native.execute();
    const functionObjects::properties& props =
        runTime.functionObjects().propsDict();
    auto nat = [&](const word& key) -> doubleScalar
    {
        return props.hasResultObjectEntry(nativeName, key)
            ? doubleScalar(props.getObjectResult<scalar>(nativeName, key))
            : std::nan("");
    };

    // --- coupledForces
    coupledForces cf(mesh, coupledForces::readDict(runTime), "coupledForcesTest");
    cf.evaluate(p, turbulence->devReff()());
    const coupledForces::coefficients& c = cf.coeffs();

    struct cmp { const char* name; doubleScalar a; doubleScalar b; bool gate; };
    const cmp rows[] =
    {
        {"Cd", c.Cd, nat("Cd"), true},
        {"Cl", c.Cl, nat("Cl"), true},
        {"Cm", c.Cm, nat("CmPitch"), true},
        {"Cd_p", c.Cd_p, nat("CdPressure"), false},
        {"Cd_v", c.Cd_v, nat("CdViscous"), false},
        {"Cl_p", c.Cl_p, nat("ClPressure"), false},
        {"Cl_v", c.Cl_v, nat("ClViscous"), false}
    };

    jsonWriter j;
    j.add("time", runTime.timeName());
    j.add("nProcs", label(UPstream::nProcs()));
    j.add("native", nativeName);
    j.add("tol", tol);
    j.add("scalarBytes", label(sizeof(scalar)));

    bool pass = true;
    doubleScalar worst = 0;
    Info<< "  coefficient  coupledForces        native               rel.diff"
        << nl;
    for (const cmp& r : rows)
    {
        // GUARD: |native| floored (D3: literal double guard)
        const doubleScalar rel =
            std::abs(r.a - r.b)/std::max(std::abs(r.b), 1e-300);
        Info<< "  " << r.name << "  " << r.a << "  " << r.b << "  " << rel
            << (r.gate ? "" : "  (reported)") << nl;
        j.add(word(r.name), r.a);
        j.add(word(r.name) + "_native", r.b);
        j.add(word(r.name) + "_relDiff", rel);
        if (r.gate)
        {
            const bool ok = std::isfinite(rel) && rel < tol;
            pass = pass && ok;
            worst = std::max(worst, std::isfinite(rel) ? rel : 1.0);
        }
    }
    j.add("maxRelDiff", worst);
    j.add("pass", pass);
    Info<< nl << "Test-coupledForces: max rel. difference (Cd, Cl, Cm) "
        << worst << (pass ? "  PASS" : "  FAIL") << nl << endl;

    if (args.found("json") && UPstream::master())
    {
        j.write(args.get<fileName>("json"));
    }

    Info<< "End" << nl << endl;
    return pass ? 0 : 1;
}


// ************************************************************************* //
