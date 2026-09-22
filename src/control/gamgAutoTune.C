/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "gamgAutoTune.H"
#include "coupledDefaults.H"
#include <algorithm>
#include <cmath>
#include <cstdio>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gamgAutoTune::gamgAutoTune
(
    const blockGAMG& gamg,
    const dictionary& blockGAMGDict
)
:
    gamg_(gamg),
    autoTune_(coupledDefaults::autoTune),
    tuneInterval_(coupledDefaults::tuneInterval),
    nPostSweepsMax_(coupledDefaults::nPostSweepsMax),
    rhoHigh_(coupledDefaults::tuneRhoHigh),
    rhoLow_(coupledDefaults::tuneRhoLow),
    rhoDemote_(coupledDefaults::tuneRhoDemote),
    rhoFail_(coupledDefaults::tuneRhoFail),
    consecutiveWindows_(coupledDefaults::tuneConsecutiveWindows),
    nPostSweepsMin_(coupledDefaults::minPostSweeps),
    userCycle_(gamg.cycleType()),
    window_(),
    nInWindow_(0),
    prevCondition_(condition::none),
    nSameCondition_(0),
    nFailWindows_(0),
    nPromoted_(0),
    promotionRefused_(false),
    kSaturated_(false),
    failed_(false),
    changed_(false),
    events_()
{
    const dictionary& d = blockGAMGDict;

    autoTune_ = d.getOrDefault<bool>("autoTune", coupledDefaults::autoTune);
    tuneInterval_ =
        d.getOrDefault<label>("tuneInterval", coupledDefaults::tuneInterval);
    nPostSweepsMax_ =
        d.getOrDefault<label>
        (
            "nPostSweepsMax",
            coupledDefaults::nPostSweepsMax
        );

    rhoHigh_ = d.getOrDefault<doubleScalar>("tuneRhoHigh", rhoHigh_);
    rhoLow_ = d.getOrDefault<doubleScalar>("tuneRhoLow", rhoLow_);
    rhoDemote_ = d.getOrDefault<doubleScalar>("tuneRhoDemote", rhoDemote_);
    rhoFail_ = d.getOrDefault<doubleScalar>("tuneRhoFail", rhoFail_);
    consecutiveWindows_ =
        d.getOrDefault<label>("tuneConsecutiveWindows", consecutiveWindows_);
    nPostSweepsMin_ = d.getOrDefault<label>("nPostSweepsMin", nPostSweepsMin_);

    // Negated comparisons also reject non-finite input
    if
    (
        !(rhoDemote_ > 0 && rhoDemote_ <= rhoLow_ && rhoLow_ < rhoHigh_
       && rhoHigh_ <= rhoFail_ && rhoFail_ < 1)
    )
    {
        FatalIOErrorInFunction(d)
            << "autoTune thresholds must satisfy 0 < tuneRhoDemote <="
            << " tuneRhoLow < tuneRhoHigh <= tuneRhoFail < 1, got "
            << rhoDemote_ << ", " << rhoLow_ << ", " << rhoHigh_ << ", "
            << rhoFail_ << exit(FatalIOError);
    }
    if (consecutiveWindows_ < 1 || nPostSweepsMin_ < 1)
    {
        FatalIOErrorInFunction(d)
            << "tuneConsecutiveWindows (" << consecutiveWindows_
            << ") and nPostSweepsMin (" << nPostSweepsMin_
            << ") must be >= 1" << exit(FatalIOError);
    }

    if (tuneInterval_ < 1)
    {
        FatalIOErrorInFunction(d)
            << "tuneInterval " << tuneInterval_ << " must be >= 1"
            << exit(FatalIOError);
    }
    if (nPostSweepsMax_ < nPostSweepsMin_)
    {
        FatalIOErrorInFunction(d)
            << "nPostSweepsMax " << nPostSweepsMax_ << " must be >= "
            << nPostSweepsMin_ << exit(FatalIOError);
    }

    window_.reserve(tuneInterval_);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::doubleScalar Foam::gamgAutoTune::median(const UList<doubleScalar>& v)
{
    List<doubleScalar> s(v);
    std::sort(s.begin(), s.end());

    const label n = s.size();
    const label mid = n/2;
    if (n % 2)
    {
        return s[mid];
    }
    return (s[mid - 1] + s[mid])/2;
}


void Foam::gamgAutoTune::addEvent
(
    const label iter,
    const std::string& what,
    const doubleScalar rhoMed
)
{
    Info<< "GAMG-tune: " << what.c_str() << " rho_med=" << rhoMed << endl;
    events_.append(tuneEvent{iter, what, rhoMed});
}


void Foam::gamgAutoTune::applyHigh(const label iter, const doubleScalar rhoMed)
{
    const label nPost = gamg_.nPostSweeps();

    if (nPost < nPostSweepsMax_)
    {
        gamg_.setNPostSweeps(nPost + 1);
        changed_ = true;
        addEvent
        (
            iter,
            "nPostSweeps " + std::to_string(nPost)
          + "->" + std::to_string(nPost + 1),
            rhoMed
        );
        return;
    }

    const blockGAMG::cycleKind c = gamg_.cycleType();

    if (c == blockGAMG::cycleKind::V || c == blockGAMG::cycleKind::F)
    {
        const blockGAMG::cycleKind next =
        (
            c == blockGAMG::cycleKind::V
          ? blockGAMG::cycleKind::F
          : blockGAMG::cycleKind::W
        );

        if (gamg_.setCycleType(next))
        {
            ++nPromoted_;
            changed_ = true;
            addEvent
            (
                iter,
                "cycleType " + blockGAMG::cycleName(c)
              + "->" + blockGAMG::cycleName(next),
                rhoMed
            );
        }
        else if (!promotionRefused_)
        {
            promotionRefused_ = true;
            addEvent(iter, "promotion-refused", rhoMed);
        }
    }
    else if (c == blockGAMG::cycleKind::K && !kSaturated_)
    {
        kSaturated_ = true;
        addEvent(iter, "kcycle-saturated", rhoMed);
    }
}


Foam::blockGAMG::cycleKind Foam::gamgAutoTune::demotionTarget() const
{
    const blockGAMG::cycleKind c = gamg_.cycleType();
    blockGAMG::cycleKind prev = c;
    if (c == blockGAMG::cycleKind::W)
    {
        prev = blockGAMG::cycleKind::F;
    }
    else if (c == blockGAMG::cycleKind::F)
    {
        prev = blockGAMG::cycleKind::V;
    }

    // Only a controller promotion is undone, never below the user's cycle
    // type (V < F < W); K is never promoted to
    if
    (
        nPromoted_ < 1
     || prev == c
     || static_cast<int>(prev) < static_cast<int>(userCycle_)
    )
    {
        return c;
    }
    return prev;
}


void Foam::gamgAutoTune::demote(const label iter, const doubleScalar rhoMed)
{
    const blockGAMG::cycleKind c = gamg_.cycleType();
    const blockGAMG::cycleKind prev = demotionTarget();

    if (prev != c && gamg_.setCycleType(prev))
    {
        --nPromoted_;
        changed_ = true;
        addEvent
        (
            iter,
            "cycleType " + blockGAMG::cycleName(c)
          + "->" + blockGAMG::cycleName(prev),
            rhoMed
        );
    }
}


void Foam::gamgAutoTune::applyLow(const label iter, const doubleScalar rhoMed)
{
    // 6.3.5: "rho < 0.3 and nPostSweeps > 1: nPostSweeps -= 1. If
    // nPostSweeps == 1 and rho < 0.2 and cycleType was promoted by the
    // controller: demote one step." The demotion is checked after the
    // decrement (reached on the 2 -> 1 step) and also on its own at
    // nPostSweeps == 1 (the low condition of evaluate() includes that case),
    // so it stays reachable once the sweeps are at the minimum.
    const label nPost = gamg_.nPostSweeps();

    if (nPost > nPostSweepsMin_)
    {
        gamg_.setNPostSweeps(nPost - 1);
        changed_ = true;
        addEvent
        (
            iter,
            "nPostSweeps " + std::to_string(nPost)
          + "->" + std::to_string(nPost - 1),
            rhoMed
        );

        // Saturation latches refer to the maximum sweep count
        promotionRefused_ = false;
        kSaturated_ = false;
    }

    if
    (
        gamg_.nPostSweeps() == nPostSweepsMin_
     && rhoMed < rhoDemote_
    )
    {
        demote(iter, rhoMed);
    }
}


void Foam::gamgAutoTune::evaluate(const label iter, const doubleScalar rhoMed)
{
    // Failure condition (checked on the state the window was run with)
    const blockGAMG::cycleKind c = gamg_.cycleType();
    if
    (
        rhoMed > rhoFail_
     && gamg_.nPostSweeps() >= nPostSweepsMax_
     && (c == blockGAMG::cycleKind::W || c == blockGAMG::cycleKind::K)
    )
    {
        ++nFailWindows_;
    }
    else
    {
        nFailWindows_ = 0;
    }

    if (nFailWindows_ >= consecutiveWindows_ && !failed_)
    {
        failed_ = true;
        WarningInFunction
            << "blockGAMG cycle efficiency failure: rho_med " << rhoMed
            << " > " << rhoFail_
            << " in " << consecutiveWindows_
            << " consecutive windows at nPostSweepsMax "
            << nPostSweepsMax_ << " with cycleType "
            << blockGAMG::cycleName(c) << endl;
        addEvent(iter, "failure", rhoMed);
    }

    // Classify the window. low: a sweep can be removed, or the sweeps are
    // at the minimum and a controller promotion can be undone (6.3.5)
    const label nPost = gamg_.nPostSweeps();
    condition cond = condition::none;
    if (rhoMed > rhoHigh_)
    {
        cond = condition::high;
    }
    else if
    (
        rhoMed < rhoLow_
     && (
            nPost > nPostSweepsMin_
         || (
                rhoMed < rhoDemote_
             && demotionTarget() != gamg_.cycleType()
            )
        )
    )
    {
        cond = condition::low;
    }

    // Hysteresis: the same condition in tuneConsecutiveWindows consecutive
    // windows; after a change the count restarts
    if (cond == prevCondition_)
    {
        ++nSameCondition_;
    }
    else
    {
        prevCondition_ = cond;
        nSameCondition_ = 1;
    }

    if
    (
        cond != condition::none
     && nSameCondition_ >= consecutiveWindows_
    )
    {
        if (cond == condition::high)
        {
            applyHigh(iter, rhoMed);
        }
        else
        {
            applyLow(iter, rhoMed);
        }
        prevCondition_ = condition::none;
        nSameCondition_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::gamgAutoTune::record(const label iter, const doubleScalar rho)
{
    changed_ = false;

    if (!autoTune_)
    {
        return false;
    }

    if (rho >= 0 && std::isfinite(rho))
    {
        window_.append(rho);
    }
    ++nInWindow_;

    if (nInWindow_ < tuneInterval_)
    {
        return false;
    }

    if (window_.size())
    {
        evaluate(iter, median(window_));
    }

    window_.clear();
    nInWindow_ = 0;

    return changed_;
}


void Foam::gamgAutoTune::writeState(dictionary& dict) const
{
    dictionary d;
    d.set("nPostSweeps", gamg_.nPostSweeps());
    d.set("cycleType", blockGAMG::cycleName(gamg_.cycleType()));
    d.set("window", List<doubleScalar>(window_));
    d.set("nInWindow", nInWindow_);
    d.set("prevCondition", label(prevCondition_));
    d.set("nSameCondition", nSameCondition_);
    d.set("nFailWindows", nFailWindows_);
    d.set("nPromoted", nPromoted_);
    d.set("promotionRefused", promotionRefused_);
    d.set("kSaturated", kSaturated_);
    d.set("failed", failed_);
    dict.set("gamgAutoTune", d);
}


void Foam::gamgAutoTune::readState(const dictionary& dict)
{
    const dictionary* dp = dict.findDict("gamgAutoTune");
    if (!dp)
    {
        Info<< "gamgAutoTune: no state in the restart dictionary, starting"
            << " from the settings" << endl;
        return;
    }
    const dictionary& d = *dp;

    // Re-apply the controller's sweeps and cycle type to the GAMG
    const label nPost = d.get<label>("nPostSweeps");
    if (nPost < nPostSweepsMin_)
    {
        FatalIOErrorInFunction(d)
            << "nPostSweeps " << nPost << " < "
            << nPostSweepsMin_ << exit(FatalIOError);
    }
    gamg_.setNPostSweeps(nPost);

    const blockGAMG::cycleKind c =
        blockGAMG::cycleFromWord(d.get<word>("cycleType"));
    if (c != gamg_.cycleType() && !gamg_.setCycleType(c))
    {
        WarningInFunction
            << "restart: cycleType " << blockGAMG::cycleName(c)
            << " refused by the coarsening-ratio rule, keeping "
            << blockGAMG::cycleName(gamg_.cycleType()) << endl;
    }

    const List<doubleScalar> w(d.get<List<doubleScalar>>("window"));
    window_.clear();
    window_.push_back(w);
    nInWindow_ = d.get<label>("nInWindow");

    const label pc = d.get<label>("prevCondition");
    if (pc < label(condition::none) || pc > label(condition::low))
    {
        FatalIOErrorInFunction(d)
            << "prevCondition " << pc << " out of range" << exit(FatalIOError);
    }
    prevCondition_ = condition(pc);
    nSameCondition_ = d.get<label>("nSameCondition");
    nFailWindows_ = d.get<label>("nFailWindows");
    nPromoted_ = d.get<label>("nPromoted");
    promotionRefused_ = d.get<bool>("promotionRefused");
    kSaturated_ = d.get<bool>("kSaturated");
    failed_ = d.get<bool>("failed");

    Info<< "gamgAutoTune: restart with nPostSweeps " << gamg_.nPostSweeps()
        << ", cycleType " << blockGAMG::cycleName(gamg_.cycleType())
        << ", " << window_.size() << " rho samples in the window" << endl;
}


std::string Foam::gamgAutoTune::eventsJson() const
{
    std::string s("[");
    forAll(events_, i)
    {
        const tuneEvent& e = events_[i];
        if (i)
        {
            s += ",";
        }

        std::string what;
        for (const char ch : e.what)
        {
            if (ch == '"' || ch == '\\')
            {
                what += '\\';
            }
            what += ch;
        }

        std::string num("null");
        const doubleScalar v = e.rhoMed;
        if (std::isfinite(v))
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.17g", v);
            num = buf;
        }

        s += "{\"iter\":" + std::to_string(e.iter)
          + ",\"what\":\"" + what
          + "\",\"rhoMed\":" + num + "}";
    }
    s += "]";
    return s;
}


Foam::dictionary Foam::gamgAutoTune::settings() const
{
    dictionary d;
    d.add("autoTune", autoTune_);
    d.add("tuneInterval", tuneInterval_);
    d.add("nPostSweepsMax", nPostSweepsMax_);
    d.add("nPostSweepsMin", nPostSweepsMin_);
    d.add("tuneRhoHigh", rhoHigh_);
    d.add("tuneRhoLow", rhoLow_);
    d.add("tuneRhoDemote", rhoDemote_);
    d.add("tuneRhoFail", rhoFail_);
    d.add("tuneConsecutiveWindows", consecutiveWindows_);
    d.add("userCycleType", blockGAMG::cycleName(userCycle_));
    return d;
}


// ************************************************************************* //
