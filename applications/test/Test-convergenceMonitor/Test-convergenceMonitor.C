/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-convergenceMonitor

Description
    Unit test of the force-coefficient window statistics and the optional
    stationary-mean (drift) stop rule of convergenceMonitor (D-045 b) on
    synthetic sequences:

    1. constant Cd, Cl, Cm: mean exact, RMS 0 and drift 0 exactly; the
       min/max window rule converges
    2. sinusoid with an integer number of periods per half-window: RMS =
       amplitude/sqrt(2), drift ~0; the min/max rule fails; the drift rule
       converges only when enabled (forceCoeffsDriftTol > 0), not with the
       default (off)
    3. linear ramp: the drift rule does not converge (drift = slope*W/2)
    4. restart round trip (writeState/readState): identical statistics
    5. the drift rule needs a full window: no convergence before

Usage
    Test-convergenceMonitor [-json <file>]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "convergenceMonitor.H"
#include "jsonWriter.H"
#include "mathematicalConstants.H"
#include <cmath>

using namespace Foam;

// Test parameters (not solver parameters)
static constexpr label window = 100;
static constexpr label period = 10;           // samples per sinusoid period
static constexpr doubleScalar amp = 0.05;     // sinusoid amplitude
static constexpr doubleScalar slope = 1e-3;   // ramp per sample
static constexpr doubleScalar driftTol = 0.01;
static constexpr doubleScalar tolExact = 1e-12;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static dictionary coupledDict(const doubleScalar dTol)
{
    dictionary c;
    c.add("mode", word("any"));
    c.add("forceCoeffsWindow", window);
    c.add("residualTol", doubleScalar(1e-12));
    if (dTol > 0)
    {
        c.add("forceCoeffsDriftTol", dTol);
    }
    dictionary d;
    d.add("convergence", c);
    return d;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("json", "file", "Write results as JSON");
    argList args(argc, argv);

    const doubleScalar R = 1;   // residual criterion never holds
    bool pass = true;
    auto check = [&](const bool ok, const char* what)
    {
        Info<< (ok ? "  ok   " : "  FAIL ") << what << endl;
        pass = pass && ok;
    };

    // 1. constant
    convergenceMonitor c1(coupledDict(driftTol));
    for (label k = 0; k < 2*window; ++k)
    {
        c1.addSample(0.3, -0.7, 0.01);
    }
    Info<< "constant sequence" << endl;
    check(c1.nCoeffs() == 3, "Cd, Cl, Cm monitored");
    bool exact = true;
    for (label i = 0; i < c1.nCoeffs(); ++i)
    {
        const convergenceMonitor::coeffStats st = c1.stats(i);
        exact = exact && st.rms == 0 && st.drift == 0 && st.n == window
             && std::abs(st.mean - c1.last(i)) <= tolExact;
    }
    check(exact, "RMS 0 and drift 0 exactly, mean = value");
    check(c1.converged(R), "converged (min/max window rule)");
    check(c1.driftConverged(), "drift rule holds");

    // 2. sinusoid, drift rule off (default) and on
    const doubleScalar w = constant::mathematical::twoPi/scalar(period);
    convergenceMonitor off(coupledDict(0));
    convergenceMonitor on(coupledDict(driftTol));
    for (label k = 0; k < 3*window; ++k)
    {
        const doubleScalar s = amp*std::sin(w*scalar(k));
        off.addSample(1 + s, 0.5 + s);
        on.addSample(1 + s, 0.5 + s);
    }
    const convergenceMonitor::coeffStats sd = on.stats(0);
    Info<< "sinusoid: mean " << sd.mean << " rms " << sd.rms
        << " drift " << sd.drift << endl;
    check(std::abs(sd.rms - amp/std::sqrt(2.0)) < 1e-9, "RMS = amp/sqrt(2)");
    check(sd.drift < 1e-9, "drift ~ 0");
    check(!off.converged(R), "drift rule off (default): not converged");
    check(on.converged(R), "drift rule on: converged");

    // 3. ramp
    convergenceMonitor ramp(coupledDict(driftTol));
    for (label k = 0; k < 2*window; ++k)
    {
        ramp.addSample(1 + slope*scalar(k), 0.5);
    }
    const convergenceMonitor::coeffStats sr = ramp.stats(0);
    Info<< "ramp: drift " << sr.drift << endl;
    check
    (
        std::abs(sr.drift - slope*scalar(window/2)) < 1e-9,
        "drift = slope*W/2"
    );
    check(!ramp.converged(R), "ramp not converged");

    // 4. restart round trip
    dictionary st;
    on.writeState(st);
    convergenceMonitor re(coupledDict(driftTol));
    re.readState(st);
    re.addSample(1.01, 0.51);
    on.addSample(1.01, 0.51);
    bool same = true;
    for (label i = 0; i < 2; ++i)
    {
        const convergenceMonitor::coeffStats a = on.stats(i);
        const convergenceMonitor::coeffStats b = re.stats(i);
        same = same && a.mean == b.mean && a.rms == b.rms
            && a.drift == b.drift && a.n == b.n;
        Info<< "  restart " << i << ": n " << a.n << " " << b.n
            << " mean " << a.mean - b.mean << " rms " << a.rms - b.rms
            << " drift " << a.drift - b.drift << endl;
    }
    check(same, "restart: identical window statistics");

    // 5. window not full
    convergenceMonitor part(coupledDict(driftTol));
    for (label k = 0; k < window - 1; ++k)
    {
        part.addSample(0.3, -0.7);
    }
    check(!part.driftConverged(), "drift rule needs a full window");
    check(part.stats(0).n == window - 1, "statistics over the samples so far");

    Info<< (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-convergenceMonitor");
        j.add("sinusoidRms", sd.rms);
        j.add("sinusoidDrift", sd.drift);
        j.add("rampDrift", sr.drift);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
