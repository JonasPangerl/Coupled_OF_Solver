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
    userCycle_(gamg.cycleType()),
    window_(),
    nInWindow_(0),
    prevCondition_(condition::none),
    nFailWindows_(0),
    nPromoted_(0),
    promotionRefused_(false),
    kSaturated_(false),
    failed_(false),
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

    if (tuneInterval_ < 1)
    {
        FatalIOErrorInFunction(d)
            << "tuneInterval " << tuneInterval_ << " must be >= 1"
            << exit(FatalIOError);
    }
    if (nPostSweepsMax_ < 1)
    {
        FatalIOErrorInFunction(d)
            << "nPostSweepsMax " << nPostSweepsMax_ << " must be >= 1"
            << exit(FatalIOError);
    }

    window_.reserve(tuneInterval_);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::gamgAutoTune::median(const UList<scalar>& v)
{
    List<scalar> s(v);
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
    const scalar rhoMed
)
{
    Info<< "GAMG-tune: " << what.c_str() << " rho_med=" << rhoMed << endl;
    events_.append(tuneEvent{iter, what, rhoMed});
}


void Foam::gamgAutoTune::applyHigh(const label iter, const scalar rhoMed)
{
    const label nPost = gamg_.nPostSweeps();

    if (nPost < nPostSweepsMax_)
    {
        gamg_.setNPostSweeps(nPost + 1);
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


void Foam::gamgAutoTune::applyLow(const label iter, const scalar rhoMed)
{
    const label nPost = gamg_.nPostSweeps();

    gamg_.setNPostSweeps(nPost - 1);
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

    if
    (
        nPost - 1 == 1
     && rhoMed < coupledDefaults::tuneRhoDemote
     && nPromoted_ > 0
    )
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

        // Never below the user's cycle type (V < F < W)
        if
        (
            prev != c
         && static_cast<int>(prev) >= static_cast<int>(userCycle_)
         && gamg_.setCycleType(prev)
        )
        {
            --nPromoted_;
            addEvent
            (
                iter,
                "cycleType " + blockGAMG::cycleName(c)
              + "->" + blockGAMG::cycleName(prev),
                rhoMed
            );
        }
    }
}


void Foam::gamgAutoTune::evaluate(const label iter, const scalar rhoMed)
{
    // Failure condition (checked on the state the window was run with)
    const blockGAMG::cycleKind c = gamg_.cycleType();
    if
    (
        rhoMed > coupledDefaults::tuneRhoFail
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

    if (nFailWindows_ >= 2 && !failed_)
    {
        failed_ = true;
        WarningInFunction
            << "blockGAMG cycle efficiency failure: rho_med " << rhoMed
            << " > " << coupledDefaults::tuneRhoFail
            << " in two consecutive windows at nPostSweepsMax "
            << nPostSweepsMax_ << " with cycleType "
            << blockGAMG::cycleName(c) << endl;
        addEvent(iter, "failure", rhoMed);
    }

    // Classify the window
    condition cond = condition::none;
    if (rhoMed > coupledDefaults::tuneRhoHigh)
    {
        cond = condition::high;
    }
    else if (rhoMed < coupledDefaults::tuneRhoLow && gamg_.nPostSweeps() > 1)
    {
        cond = condition::low;
    }

    // Hysteresis: the same condition in two consecutive windows
    if (cond != condition::none && cond == prevCondition_)
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
    }
    else
    {
        prevCondition_ = cond;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gamgAutoTune::record(const label iter, const scalar rho)
{
    if (!autoTune_)
    {
        return;
    }

    if (rho >= 0 && std::isfinite(rho))
    {
        window_.append(rho);
    }
    ++nInWindow_;

    if (nInWindow_ < tuneInterval_)
    {
        return;
    }

    if (window_.size())
    {
        evaluate(iter, median(window_));
    }

    window_.clear();
    nInWindow_ = 0;
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
    d.add("userCycleType", blockGAMG::cycleName(userCycle_));
    return d;
}


// ************************************************************************* //
