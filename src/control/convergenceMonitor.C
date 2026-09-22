/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "convergenceMonitor.H"
#include "coupledDefaults.H"
#include "functionObjectList.H"
#include "functionObjectProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceMonitor::convergenceMonitor(const dictionary& coupledDict)
:
    mode_("any"),
    window_(coupledDefaults::forceCoeffsWindow),
    forceTol_(coupledDefaults::forceCoeffsTol),
    residualTol_(coupledDefaults::residualTol),
    Cd_(),
    Cl_(),
    foName_()
{
    const dictionary& d = coupledDict.subOrEmptyDict("convergence");
    mode_ = d.getOrDefault<word>("mode", "any");
    window_ = d.getOrDefault<label>("forceCoeffsWindow", window_);
    forceTol_ = d.getOrDefault<scalar>("forceCoeffsTol", forceTol_);
    residualTol_ = d.getOrDefault<scalar>("residualTol", residualTol_);

    if (mode_ != "any" && mode_ != "all")
    {
        FatalIOErrorInFunction(d)
            << "convergence.mode must be any or all, not " << mode_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::convergenceMonitor::windowConverged
(
    const DynamicList<scalar>& h
) const
{
    if (window_ <= 0 || h.size() < window_)
    {
        return false;
    }

    scalar mn = GREAT, mx = -GREAT, sum = 0;
    for (label i = h.size() - window_; i < h.size(); ++i)
    {
        mn = min(mn, h[i]);
        mx = max(mx, h[i]);
        sum += h[i];
    }
    const scalar mean = sum/scalar(window_);

    return (mx - mn) <= forceTol_*mag(mean);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
                props.hasResultObjectEntry(name, "Cd")
             && props.hasResultObjectEntry(name, "Cl")
            )
            {
                foName_ = name;
                Info<< "convergenceMonitor: using force coefficients of "
                    << foName_ << endl;
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
    const scalar Cd = props.getObjectResult<scalar>(foName_, "Cd");
    const scalar Cl = props.getObjectResult<scalar>(foName_, "Cl");

    if (Cd_.size() && Cd == Cd_.last() && Cl == Cl_.last())
    {
        return;
    }

    Cd_.append(Cd);
    Cl_.append(Cl);
}


void Foam::convergenceMonitor::writeState(dictionary& dict) const
{
    // The last window samples are all the criterion needs
    const label n = (window_ > 0 ? min(window_, Cd_.size()) : 0);
    const label start = Cd_.size() - n;

    dictionary d;
    d.set("Cd", scalarList(SubList<scalar>(Cd_, n, start)));
    d.set("Cl", scalarList(SubList<scalar>(Cl_, n, start)));
    dict.set("convergenceMonitor", d);
}


void Foam::convergenceMonitor::readState(const dictionary& dict)
{
    const dictionary* dp = dict.findDict("convergenceMonitor");
    if (!dp)
    {
        return;
    }

    const scalarList Cd(dp->get<scalarList>("Cd"));
    const scalarList Cl(dp->get<scalarList>("Cl"));
    if (Cd.size() != Cl.size())
    {
        FatalIOErrorInFunction(*dp)
            << "Cd and Cl histories differ in size" << exit(FatalIOError);
    }

    Cd_.clear();
    Cl_.clear();
    Cd_.push_back(Cd);
    Cl_.push_back(Cl);
}


bool Foam::convergenceMonitor::converged(const scalar R) const
{
    const bool resOk = R < residualTol_;

    if (!haveForces())
    {
        return resOk;
    }

    const bool forceOk = windowConverged(Cd_) && windowConverged(Cl_);

    if (mode_ == "all")
    {
        return resOk && forceOk;
    }
    return resOk || forceOk;
}


Foam::scalar Foam::convergenceMonitor::forceCriterionRatio() const
{
    if (foName_.empty() || window_ <= 0 || forceTol_ <= 0)
    {
        return -1;
    }
    scalar worst = 0;
    for (const DynamicList<scalar>* hp : {&Cd_, &Cl_})
    {
        const DynamicList<scalar>& h = *hp;
        if (h.size() < window_)
        {
            return -1;
        }
        scalar mn = GREAT, mx = -GREAT, sum = 0;
        for (label i = h.size() - window_; i < h.size(); ++i)
        {
            mn = min(mn, h[i]);
            mx = max(mx, h[i]);
            sum += h[i];
        }
        // GUARD: |mean| floored at VSMALL
        const scalar den = forceTol_*max(mag(sum/scalar(window_)), VSMALL);
        worst = max(worst, (mx - mn)/den);
    }
    return worst;
}


void Foam::convergenceMonitor::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("mode", mode_);
    d.add("forceCoeffsWindow", window_);
    d.add("forceCoeffsTol", forceTol_);
    d.add("residualTol", residualTol_);
    dict.add("convergence", d);
}


// ************************************************************************* //
