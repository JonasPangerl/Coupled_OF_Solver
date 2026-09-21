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

    if (!foName_.empty())
    {
        Cd_.append(props.getObjectResult<scalar>(foName_, "Cd"));
        Cl_.append(props.getObjectResult<scalar>(foName_, "Cl"));
    }
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
