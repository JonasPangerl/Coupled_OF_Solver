/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "startupControl.H"
#include "coupledDefaults.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::startupControl::startupControl
(
    const dictionary& coupledDict
)
:
    mode_(mode::hybrid),
    modeName_(coupledDict.getOrDefault<word>
    (
        "startupMode",
        word(coupledDefaults::startupMode)
    )),
    upwindIters_(coupledDict.getOrDefault<label>
    (
        "startupUpwindIters",
        coupledDefaults::startupUpwindIters
    )),
    switchR_(coupledDict.getOrDefault<doubleScalar>
    (
        "startupSwitchR",
        coupledDefaults::startupSwitchR
    )),
    rampStart_(coupledDict.getOrDefault<label>
    (
        "startupRampStart",
        coupledDefaults::startupRampStart
    )),
    rampLength_(coupledDict.getOrDefault<label>
    (
        "startupRampLength",
        coupledDefaults::startupRampLength
    )),
    rampStartMax_(coupledDict.getOrDefault<label>
    (
        "startupRampStartMax",
        coupledDefaults::startupRampStartMax
    )),
    stagnationTrigger_(coupledDict.getOrDefault<bool>
    (
        "startupStagnationTrigger",
        coupledDefaults::startupStagnationTrigger
    )),
    stagnationFactor_(coupledDict.getOrDefault<doubleScalar>
    (
        "startupStagnationFactor",
        coupledDefaults::startupStagnationFactor
    )),
    stagnationWindow_(coupledDict.getOrDefault<label>
    (
        "startupStagnationWindow",
        coupledDefaults::startupStagnationWindow
    )),
    fastFactor_(coupledDict.getOrDefault<doubleScalar>
    (
        "startupFastFactor",
        coupledDefaults::startupFastFactor
    )),
    developedTol_(coupledDict.getOrDefault<doubleScalar>
    (
        "startupDevelopedTol",
        coupledDefaults::startupDevelopedTol
    )),
    rampStartIter_(-1),
    trigger_(),
    Rhist_(),
    probing_(false),
    full_(false)
{
    if (modeName_ == "hybrid")
    {
        mode_ = mode::hybrid;
    }
    else if (modeName_ == "upwind")
    {
        mode_ = mode::upwind;
    }
    else if (modeName_ == "none")
    {
        mode_ = mode::none;
    }
    else
    {
        FatalIOErrorInFunction(coupledDict)
            << "Unknown startupMode " << modeName_
            << ", valid: hybrid upwind none" << exit(FatalIOError);
    }

    if
    (
        rampStart_ < 1 || rampLength_ < 0 || rampStartMax_ < rampStart_
     || stagnationWindow_ < 1
     || !(stagnationFactor_ > 0) || !(fastFactor_ > 0)
     || !(developedTol_ >= 0)
    )
    {
        FatalIOErrorInFunction(coupledDict)
            << "Invalid start-up settings: require startupRampStart >= 1,"
            << " startupRampLength >= 0, startupRampStartMax >="
            << " startupRampStart, startupStagnationWindow >= 1,"
            << " startupStagnationFactor > 0, startupFastFactor > 0,"
            << " startupDevelopedTol >= 0" << exit(FatalIOError);
    }

    if (mode_ == mode::none)
    {
        full_ = true;
        rampStartIter_ = 0;
        trigger_ = "none";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::startupControl::beta(const label m) const
{
    if (full_ || (probing_ && m == 1))
    {
        return 1;
    }
    if (rampStartIter_ < 0)
    {
        return 0;
    }
    if (mode_ == mode::upwind || rampLength_ == 0)
    {
        return (m > rampStartIter_ ? 1 : 0);
    }
    // GUARD: rampLength_ > 0 here
    const scalar b = scalar(m - rampStartIter_)/scalar(rampLength_);
    return min(scalar(1), max(scalar(0), b));
}


Foam::label Foam::startupControl::rampEndIter() const
{
    if (full_)
    {
        return (mode_ == mode::none || trigger_ == "developedStart") ? 1 : -1;
    }
    if (rampStartIter_ < 0)
    {
        return -1;
    }
    if (mode_ == mode::upwind || rampLength_ == 0)
    {
        return rampStartIter_ + 1;
    }
    return rampStartIter_ + rampLength_;
}


bool Foam::startupControl::startProbe(const bool candidate)
{
    probing_ = (mode_ != mode::none && candidate && !full_);
    return probing_;
}


bool Foam::startupControl::decideDeveloped(const doubleScalar rU, const doubleScalar rp)
{
    probing_ = false;
    if (max(rU, rp) < developedTol_)
    {
        full_ = true;
        rampStartIter_ = 0;
        trigger_ = "developedStart";
        return true;
    }
    return false;
}


bool Foam::startupControl::update(const label n, const doubleScalar R)
{
    if (full_ || rampStartIter_ >= 0)
    {
        return false;
    }

    // R of the last W+1 accepted iterations, oldest first
    const label W = stagnationWindow_;
    Rhist_.append(R);
    if (Rhist_.size() > W + 1)
    {
        Rhist_ = List<doubleScalar>(SubList<doubleScalar>(Rhist_, W + 1, Rhist_.size() - W - 1));
    }
    const bool haveWindow = (Rhist_.size() == W + 1);
    const doubleScalar Rold = Rhist_.first();

    word trig;

    if (R < switchR_)
    {
        trig = "relativeR";
    }
    else if (mode_ == mode::upwind)
    {
        if (n >= upwindIters_)
        {
            trig = "fixed";
        }
    }
    else if (n < rampStart_)
    {
        // Upwind no longer productive: less than a (1 - factor) reduction
        // over the window
        if (stagnationTrigger_ && haveWindow && R > stagnationFactor_*Rold)
        {
            trig = "stagnation";
        }
    }
    else
    {
        const bool fastDrop = haveWindow && R < fastFactor_*Rold;
        if (rampStartMax_ <= rampStart_ || (!fastDrop && n <= rampStart_))
        {
            trig = "fixed";
        }
        else if (!fastDrop)
        {
            trig = "delayedFastDrop";
        }
        else if (n >= rampStartMax_)
        {
            trig = "cap";
        }
    }

    if (trig.empty())
    {
        return false;
    }

    rampStartIter_ = n;
    trigger_ = trig;
    Rhist_.clear();
    return true;
}


void Foam::startupControl::writeState(dictionary& dict) const
{
    dict.set("startupRampStartIter", rampStartIter_);
    dict.set("startupTrigger", trigger_.empty() ? word("pending") : trigger_);
    dict.set("startupFull", full_);
    dict.set("startupRHistory", Rhist_);
}


void Foam::startupControl::readState(const dictionary& dict, const label iter)
{
    probing_ = false;
    if (dict.found("startupRampStartIter"))
    {
        rampStartIter_ = dict.get<label>("startupRampStartIter");
        trigger_ = dict.get<word>("startupTrigger");
        if (trigger_ == "pending")
        {
            trigger_.clear();
        }
        full_ = dict.getOrDefault<bool>("startupFull", false);
        Rhist_ = dict.getOrDefault<List<doubleScalar>>("startupRHistory", List<doubleScalar>());
    }
    else if (dict.getOrDefault<bool>("startupDone", false))
    {
        // State written before D-048: start-up finished
        full_ = true;
        rampStartIter_ = iter;
        trigger_ = "restartDone";
    }
}


void Foam::startupControl::writeSettings(dictionary& dict) const
{
    dict.add("startupMode", modeName_);
    dict.add("startupUpwindIters", upwindIters_);
    dict.add("startupSwitchR", switchR_);
    dict.add("startupRampStart", rampStart_);
    dict.add("startupRampLength", rampLength_);
    dict.add("startupRampStartMax", rampStartMax_);
    dict.add("startupStagnationTrigger", stagnationTrigger_);
    dict.add("startupStagnationFactor", stagnationFactor_);
    dict.add("startupStagnationWindow", stagnationWindow_);
    dict.add("startupFastFactor", fastFactor_);
    dict.add("startupDevelopedTol", developedTol_);
}


// ************************************************************************* //
