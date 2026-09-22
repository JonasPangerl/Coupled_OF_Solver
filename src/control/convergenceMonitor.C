/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "convergenceMonitor.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "precisionProfile.H"
#include "functionObjectList.H"
#include "functionObjectProperties.H"
#include <limits>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceMonitor::convergenceMonitor(const dictionary& coupledDict)
:
    mode_(coupledDefaults::convergenceMode),
    window_(coupledDefaults::forceCoeffsWindow),
    forceTol_(coupledDefaults::forceCoeffsTol),
    residualTol_(precisionProfile::current().residualTol),   // D7
    rmsWindow_(coupledDefaults::forceCoeffsWindow),
    driftTol_(coupledDefaults::forceCoeffsDriftTol),
    driftAbs_(coupledDefaults::forceCoeffsDriftAbs),
    Cd_(),
    Cl_(),
    Cm_(),
    haveCm_(false),
    foName_()
{
    const dictionary& d = coupledDict.subOrEmptyDict("convergence");
    mode_ = d.getOrDefault<word>("mode", mode_);
    // Force-coefficient function object (D-066): empty = the first one
    // that provides Cd and Cl
    foRequested_ = d.getOrDefault<word>("forceCoeffs", word::null);
    window_ = d.getOrDefault<label>("forceCoeffsWindow", window_);
    forceTol_ = d.getOrDefault<doubleScalar>("forceCoeffsTol", forceTol_);
    residualTol_ = d.getOrDefault<doubleScalar>("residualTol", residualTol_);
    // Statistics/drift window: default = forceCoeffsWindow
    rmsWindow_ = d.getOrDefault<label>("forceCoeffsRmsWindow", window_);
    driftTol_ = d.getOrDefault<doubleScalar>("forceCoeffsDriftTol", driftTol_);
    driftAbs_ = d.getOrDefault<doubleScalar>("forceCoeffsDriftAbs", driftAbs_);

    if (mode_ != "any" && mode_ != "all")
    {
        FatalIOErrorInFunction(d)
            << "convergence.mode must be any or all, not " << mode_
            << exit(FatalIOError);
    }
    if (rmsWindow_ < 1 || driftTol_ < 0 || driftAbs_ < 0)
    {
        FatalIOErrorInFunction(d)
            << "convergence: forceCoeffsRmsWindow must be >= 1 and"
            << " forceCoeffsDriftTol, forceCoeffsDriftAbs >= 0 (got "
            << rmsWindow_ << ", " << driftTol_ << ", " << driftAbs_ << ")"
            << exit(FatalIOError);
    }
    if (driftTol_ > 0 && rmsWindow_ < 2)
    {
        FatalIOErrorInFunction(d)
            << "convergence: the drift rule (forceCoeffsDriftTol > 0) needs"
            << " forceCoeffsRmsWindow >= 2, got " << rmsWindow_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::convergenceMonitor::windowConverged
(
    const DynamicList<doubleScalar>& h
) const
{
    if (window_ <= 0 || h.size() < window_)
    {
        return false;
    }

    doubleScalar mn = cfGreat<doubleScalar>(), mx = -cfGreat<doubleScalar>(), sum = 0;
    for (label i = h.size() - window_; i < h.size(); ++i)
    {
        mn = min(mn, h[i]);
        mx = max(mx, h[i]);
        sum += h[i];
    }
    const doubleScalar mean = sum/doubleScalar(window_);

    return (mx - mn) <= forceTol_*mag(mean);
}


Foam::convergenceMonitor::coeffStats Foam::convergenceMonitor::windowStats
(
    const DynamicList<doubleScalar>& h
) const
{
    coeffStats st;
    st.n = min(rmsWindow_, h.size());
    if (st.n < 1)
    {
        return st;
    }
    const label start = h.size() - st.n;

    // Shifted by the first sample of the window: exact for a constant
    // sequence (RMS and drift exactly 0) and free of cancellation
    const doubleScalar x0 = h[start];
    doubleScalar sum = 0;
    for (label i = start; i < h.size(); ++i)
    {
        sum += h[i] - x0;
    }
    const doubleScalar dm = sum/doubleScalar(st.n);
    st.mean = x0 + dm;

    doubleScalar var = 0;
    for (label i = start; i < h.size(); ++i)
    {
        var += sqr((h[i] - x0) - dm);
    }
    st.rms = std::sqrt(var/doubleScalar(st.n));

    if (st.n >= 2)
    {
        // Older half: the first n/2 samples of the window
        const label nOld = st.n/2;
        doubleScalar sOld = 0, sNew = 0;
        for (label i = start; i < start + nOld; ++i)
        {
            sOld += h[i] - x0;
        }
        for (label i = start + nOld; i < h.size(); ++i)
        {
            sNew += h[i] - x0;
        }
        st.drift = mag(sOld/doubleScalar(nOld) - sNew/doubleScalar(st.n - nOld));
    }
    return st;
}


const Foam::DynamicList<Foam::doubleScalar>& Foam::convergenceMonitor::history
(
    const label i
) const
{
    return (i == 0 ? Cd_ : (i == 1 ? Cl_ : Cm_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const char* Foam::convergenceMonitor::coeffName(const label i)
{
    return (i == 0 ? "Cd" : (i == 1 ? "Cl" : "Cm"));
}


Foam::doubleScalar Foam::convergenceMonitor::last(const label i) const
{
    const DynamicList<doubleScalar>& h = history(i);
    return (h.size() ? h.last() : 0);
}


Foam::convergenceMonitor::coeffStats Foam::convergenceMonitor::stats
(
    const label i
) const
{
    return windowStats(history(i));
}


bool Foam::convergenceMonitor::driftConverged() const
{
    if (!driftEnabled() || !haveForces())
    {
        return false;
    }
    for (label i = 0; i < nCoeffs(); ++i)
    {
        const DynamicList<doubleScalar>& h = history(i);
        if (h.size() < rmsWindow_)
        {
            return false;
        }
        const coeffStats st = windowStats(h);
        if (st.drift > max(driftTol_*mag(st.mean), driftAbs_))
        {
            return false;
        }
    }
    return true;
}


void Foam::convergenceMonitor::addSample(const doubleScalar Cd, const doubleScalar Cl)
{
    if (foName_.empty())
    {
        foName_ = "samples";
    }
    Cd_.append(Cd);
    Cl_.append(Cl);
}


void Foam::convergenceMonitor::addSample
(
    const doubleScalar Cd,
    const doubleScalar Cl,
    const doubleScalar Cm
)
{
    addSample(Cd, Cl);
    Cm_.append(Cm);
    haveCm_ = true;
}


void Foam::convergenceMonitor::record(const Time& runTime)
{
    const functionObjects::properties& props =
        runTime.functionObjects().propsDict();

    if (foName_.empty())
    {
        for (const word& name : props.objectResultNames())
        {
            if
            (
                (foRequested_.empty() || name == foRequested_)
             && props.hasResultObjectEntry(name, "Cd")
             && props.hasResultObjectEntry(name, "Cl")
            )
            {
                foName_ = name;
                haveCm_ = props.hasResultObjectEntry(name, "CmPitch");
                Info<< "convergenceMonitor: using force coefficients of "
                    << foName_ << (haveCm_ ? " (Cd, Cl, CmPitch)" : "")
                    << endl;
                break;
            }
        }
    }

    if (foName_.empty())
    {
        return;
    }

    // One sample per function-object evaluation. The result stored in the
    // properties is the one of the last execution; it is not renewed on
    // iterations where the object does not execute (executeInterval > 1),
    // and after a rolled-back or skipped iteration the object evaluates the
    // restored, already recorded state. Such repeats are recognised by an
    // unchanged (Cd, Cl) pair and not recorded again: a genuinely new
    // evaluation of a changed flow field reproduces both values bit for
    // bit only if the forces did not change at all.
    const doubleScalar Cd = props.getObjectResult<scalar>(foName_, "Cd");
    const doubleScalar Cl = props.getObjectResult<scalar>(foName_, "Cl");

    if (Cd_.size() && Cd == Cd_.last() && Cl == Cl_.last())
    {
        return;
    }

    if (haveCm_)
    {
        addSample(Cd, Cl, props.getObjectResult<scalar>(foName_, "CmPitch"));
    }
    else
    {
        addSample(Cd, Cl);
    }
}


void Foam::convergenceMonitor::record
(
    const scalar Cd,
    const scalar Cl,
    const scalar Cm,
    const word& source
)
{
    if (foName_.empty())
    {
        foName_ = source;
        haveCm_ = true;
        Info<< "convergenceMonitor: using force coefficients of " << foName_
            << " (Cd, Cl, Cm)" << endl;
    }

    if (Cd_.size() && Cd == Cd_.last() && Cl == Cl_.last())
    {
        return;
    }
    addSample(Cd, Cl, Cm);
}


void Foam::convergenceMonitor::writeState(dictionary& dict) const
{
    // The last window samples are all the criteria and the window
    // statistics need (mean/RMS continue seamlessly after a restart)
    // Round-trip precision: a dictionary entry is formatted with
    // IOstream::defaultPrecision() when it is created (6 digits by default)
    const unsigned oldPrecision = IOstream::defaultPrecision
    (
        std::numeric_limits<doubleScalar>::max_digits10
    );

    const label keep = max(window_, rmsWindow_);
    const label n = (keep > 0 ? min(keep, Cd_.size()) : 0);
    const label start = Cd_.size() - n;

    dictionary d;
    d.set("Cd", List<doubleScalar>(SubList<doubleScalar>(Cd_, n, start)));
    d.set("Cl", List<doubleScalar>(SubList<doubleScalar>(Cl_, n, start)));
    if (haveCm_)
    {
        const label nm = (keep > 0 ? min(keep, Cm_.size()) : 0);
        d.set("Cm", List<doubleScalar>(SubList<doubleScalar>(Cm_, nm, Cm_.size() - nm)));
    }
    dict.set("convergenceMonitor", d);

    IOstream::defaultPrecision(oldPrecision);
}


void Foam::convergenceMonitor::readState(const dictionary& dict)
{
    const dictionary* dp = dict.findDict("convergenceMonitor");
    if (!dp)
    {
        return;
    }

    const List<doubleScalar> Cd(dp->get<List<doubleScalar>>("Cd"));
    const List<doubleScalar> Cl(dp->get<List<doubleScalar>>("Cl"));
    if (Cd.size() != Cl.size())
    {
        FatalIOErrorInFunction(*dp)
            << "Cd and Cl histories differ in size" << exit(FatalIOError);
    }

    Cd_.clear();
    Cl_.clear();
    Cd_.push_back(Cd);
    Cl_.push_back(Cl);

    Cm_.clear();
    List<doubleScalar> Cm;
    if (dp->readIfPresent("Cm", Cm))
    {
        Cm_.push_back(Cm);
        haveCm_ = true;
    }
}


bool Foam::convergenceMonitor::converged(const doubleScalar R) const
{
    const bool resOk = R < residualTol_;

    if (!haveForces())
    {
        return resOk;
    }

    // (ii): the min/max window rule, or the drift rule if enabled
    const bool forceOk =
        (windowConverged(Cd_) && windowConverged(Cl_))
     || driftConverged();

    if (mode_ == "all")
    {
        return resOk && forceOk;
    }
    return resOk || forceOk;
}


Foam::doubleScalar Foam::convergenceMonitor::forceCriterionRatio() const
{
    if (foName_.empty() || window_ <= 0 || forceTol_ <= 0)
    {
        return -1;
    }
    doubleScalar worst = 0;
    for (const DynamicList<doubleScalar>* hp : {&Cd_, &Cl_})
    {
        const DynamicList<doubleScalar>& h = *hp;
        if (h.size() < window_)
        {
            return -1;
        }
        doubleScalar mn = cfGreat<doubleScalar>(), mx = -cfGreat<doubleScalar>(), sum = 0;
        for (label i = h.size() - window_; i < h.size(); ++i)
        {
            mn = min(mn, h[i]);
            mx = max(mx, h[i]);
            sum += h[i];
        }
        // GUARD: |mean| floored at cfVSmall
        const doubleScalar den = forceTol_*max(mag(sum/doubleScalar(window_)), cfVSmall<doubleScalar>());
        worst = max(worst, (mx - mn)/den);
    }
    return worst;
}


void Foam::convergenceMonitor::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("mode", mode_);
    d.add("forceCoeffs", foRequested_);
    d.add("forceCoeffsWindow", window_);
    d.add("forceCoeffsTol", forceTol_);
    d.add("residualTol", residualTol_);
    d.add("forceCoeffsRmsWindow", rmsWindow_);
    d.add("forceCoeffsDriftTol", driftTol_);
    d.add("forceCoeffsDriftAbs", driftAbs_);
    dict.add("convergence", d);
}


// ************************************************************************* //
