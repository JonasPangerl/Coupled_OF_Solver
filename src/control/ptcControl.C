/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "ptcControl.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include "DynamicList.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptcControl::ptcControl
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    strategy_(strategy::mRDM),
    strategyName_("mRDM"),
    CFL0_(coupledDefaults::CFL0),
    CFLmin_(coupledDefaults::CFLmin),
    CFLmax_(coupledDefaults::CFLmax),
    gamma_(coupledDefaults::ptcGamma),
    betaMax_(coupledDefaults::betaMax),
    betaExp_(coupledDefaults::betaExp),
    nHold_(coupledDefaults::nHold),
    continuityFactor_(coupledDefaults::ptcContinuityFactor),
    failCeiling_(coupledDefaults::ptcFailCeiling),
    ceilingRelax_(coupledDefaults::ptcCeilingRelax),
    ceiling_(coupledDefaults::CFLmax),
    localLimit_(coupledDefaults::localLimitEnabled),
    fLoc_(coupledDefaults::fLoc),
    localImplicit_(coupledDefaults::localImplicit),
    localMemory_(coupledDefaults::localMemory),
    localRecovery_(coupledDefaults::localRecovery),
    localHold_(coupledDefaults::localHold),
    localStickyAfter_(coupledDefaults::localStickyAfter),
    CFL_(coupledDefaults::CFL0),
    Rprev_(-1),
    holdRemaining_(0),
    pendingCFL_(-1),
    sampleCFL_(-1),
    sampleCost_(-1),
    haveSample_(false),
    effPrev_(-1),
    effPrevCFL_(-1),
    effBest_(-1),
    sampleLin_(-1),
    sampleSolve_(-1),
    effShare_(coupledDefaults::ptcEffSolveShare),
    shareSm_(-1),
    effWindow_(coupledDefaults::ptcEffWindow),
    effExponent_(coupledDefaults::ptcEffExponent),
    winN_(-1),
    winCFL_(-1),
    winShare_(0),
    winCost_(0),
    effPrevLin_(-1),
    effTol_(coupledDefaults::ptcEffTol),
    effDir_(1),
    effStep_(coupledDefaults::betaExp),
    effStreak_(0),
    reason_(),
    effLast_(-1),
    nLocalLimited_(0),
    nLocalThrottled_(0),
    nLocalSticky_(0),
    localFactor_(),
    localFactor0_(),
    localLimitedNow_(),
    localHoldRemaining_(),
    localCount_()
{
    const dictionary& d = coupledDict.subOrEmptyDict("ptc");

    strategyName_ = d.getOrDefault<word>
    (
        "cflStrategy",
        word(coupledDefaults::cflStrategy)
    );
    if (strategyName_ == "mRDM")
    {
        strategy_ = strategy::mRDM;
    }
    else if (strategyName_ == "EXP")
    {
        strategy_ = strategy::EXP;
    }
    else if (strategyName_ == "SER")
    {
        strategy_ = strategy::SER;
    }
    else if (strategyName_ == "EFF")
    {
        strategy_ = strategy::EFF;
    }
    else
    {
        FatalIOErrorInFunction(d)
            << "Unknown cflStrategy " << strategyName_
            << ", valid: mRDM EXP SER EFF" << exit(FatalIOError);
    }

    CFL0_ = d.getOrDefault<doubleScalar>("CFL0", coupledDefaults::CFL0);
    CFLmin_ = d.getOrDefault<doubleScalar>("CFLmin", coupledDefaults::CFLmin);
    CFLmax_ = d.getOrDefault<doubleScalar>("CFLmax", coupledDefaults::CFLmax);
    gamma_ = d.getOrDefault<doubleScalar>("gamma", coupledDefaults::ptcGamma);
    betaMax_ = d.getOrDefault<doubleScalar>("betaMax", coupledDefaults::betaMax);
    betaExp_ = d.getOrDefault<doubleScalar>("betaExp", coupledDefaults::betaExp);
    nHold_ = d.getOrDefault<label>("nHold", coupledDefaults::nHold);
    effStep_ = betaExp_;
    continuityFactor_ = d.getOrDefault<doubleScalar>
    (
        "continuityFactor",
        coupledDefaults::ptcContinuityFactor
    );
    // Negated comparison also rejects non-finite input
    if (!(continuityFactor_ >= 0) || !std::isfinite(continuityFactor_))
    {
        FatalIOErrorInFunction(d)
            << "ptc.continuityFactor must be a finite value >= 0 (0 = off),"
            << " got " << continuityFactor_ << exit(FatalIOError);
    }
    failCeiling_ = d.getOrDefault<doubleScalar>
    (
        "failCeiling",
        coupledDefaults::ptcFailCeiling
    );
    ceilingRelax_ = d.getOrDefault<doubleScalar>
    (
        "ceilingRelax",
        coupledDefaults::ptcCeilingRelax
    );
    if
    (
        !(failCeiling_ >= 0 && failCeiling_ < 1)
     || !(ceilingRelax_ >= 1) || !std::isfinite(ceilingRelax_)
    )
    {
        FatalIOErrorInFunction(d)
            << "ptc.failCeiling must be in [0, 1) (0 = off) and"
            << " ptc.ceilingRelax a finite value >= 1, got " << failCeiling_
            << " and " << ceilingRelax_ << exit(FatalIOError);
    }
    ceiling_ = CFLmax_;
    effTol_ = d.getOrDefault<doubleScalar>("effTol", coupledDefaults::ptcEffTol);
    effShare_ = d.getOrDefault<doubleScalar>
    (
        "effSolveShare",
        coupledDefaults::ptcEffSolveShare
    );
    effWindow_ = d.getOrDefault<label>("effWindow", coupledDefaults::ptcEffWindow);
    effExponent_ = d.getOrDefault<doubleScalar>
    (
        "effExponent",
        coupledDefaults::ptcEffExponent
    );
    if (effWindow_ < 1 || !(effExponent_ > 0))
    {
        FatalIOErrorInFunction(d)
            << "ptc.effWindow must be >= 1 and ptc.effExponent > 0, got "
            << effWindow_ << " and " << effExponent_ << exit(FatalIOError);
    }
    if (!(effShare_ > 0 && effShare_ < 1))
    {
        FatalIOErrorInFunction(d)
            << "ptc.effSolveShare must be in (0, 1), got " << effShare_
            << exit(FatalIOError);
    }
    if (!(effTol_ >= 0 && effTol_ < 1))
    {
        FatalIOErrorInFunction(d)
            << "ptc.effTol must be in [0, 1), got " << effTol_
            << exit(FatalIOError);
    }

    const dictionary& ll = coupledDict.subOrEmptyDict("localLimit");
    localLimit_ =
        ll.getOrDefault<bool>("enabled", coupledDefaults::localLimitEnabled);
    fLoc_ = ll.getOrDefault<scalar>("fLoc", coupledDefaults::fLoc);
    localImplicit_ =
        ll.getOrDefault<bool>("implicit", coupledDefaults::localImplicit);
    localMemory_ =
        ll.getOrDefault<bool>("memory", coupledDefaults::localMemory);
    localRecovery_ =
        ll.getOrDefault<scalar>("localRecovery", coupledDefaults::localRecovery);
    localHold_ = ll.getOrDefault<label>("localHold", coupledDefaults::localHold);
    localStickyAfter_ =
        ll.getOrDefault<label>
        (
            "localStickyAfter",
            coupledDefaults::localStickyAfter
        );

    // Negated comparisons also reject non-finite input
    if (!(localRecovery_ >= 1) || !std::isfinite(localRecovery_))
    {
        FatalIOErrorInFunction(ll)
            << "localLimit.localRecovery must be a finite value >= 1, got "
            << localRecovery_ << exit(FatalIOError);
    }
    if (localHold_ < 0 || localStickyAfter_ < 0)
    {
        FatalIOErrorInFunction(ll)
            << "localLimit.localHold and localStickyAfter must be >= 0, got "
            << localHold_ << " and " << localStickyAfter_
            << exit(FatalIOError);
    }

    if (localLimit_ && localMemory_)
    {
        localFactor_.resize(mesh_.nCells(), scalar(1));
        localFactor0_.resize(mesh_.nCells(), scalar(1));
        localLimitedNow_.resize(mesh_.nCells(), false);
        localHoldRemaining_.resize(mesh_.nCells(), 0);
        localCount_.resize(mesh_.nCells(), 0);
    }
    else
    {
        localMemory_ = false;
    }

    CFL_ = CFL0_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::ptcControl::rDeltaTV
(
    const surfaceScalarField& phi,
    const volScalarField& nuEff,
    const scalarField& cflFactor
) const
{
    const scalarField& V = mesh_.V();
    auto tr = tmp<scalarField>::New(mesh_.nCells(), Zero);
    scalarField& r = tr.ref();

    // 1/2 sum_f |phi_f|
    scalarField sumPhi(mesh_.nCells(), Zero);
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, facei)
    {
        const scalar a = mag(phi[facei]);
        sumPhi[own[facei]] += a;
        sumPhi[nei[facei]] += a;
    }
    forAll(phi.boundaryField(), patchi)
    {
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        const scalarField& pp = phi.boundaryField()[patchi];
        forAll(fc, pf)
        {
            sumPhi[fc[pf]] += mag(pp[pf]);
        }
    }

    const scalarField& nu = nuEff.primitiveField();

    forAll(r, celli)
    {
        // GUARD: pow argument V > 0; denominator >= cfVSmall (5.4)
        const scalar lambda = scalar
        (
            0.5*sumPhi[celli]
          + nu[celli]*std::cbrt(max(V[celli], cfVSmall<scalar>()))
        );
        const scalar cfl = max(scalar(CFL_*cflFactor[celli]), cfVSmall<scalar>());
        r[celli] = max(lambda, cfVSmall<scalar>())/cfl;
    }

    return tr;
}


void Foam::ptcControl::beginIteration()
{
    reason_.clear();

    if (!localMemory_)
    {
        return;
    }

    // Recovery of the memory factors: after the per-cell hold, and never
    // for sticky cells (limited in localStickyAfter iterations)
    forAll(localFactor_, celli)
    {
        if (localHoldRemaining_[celli] > 0)
        {
            --localHoldRemaining_[celli];
        }
        else if
        (
            localFactor_[celli] < 1
         && !(localStickyAfter_ > 0 && localCount_[celli] >= localStickyAfter_)
        )
        {
            localFactor_[celli] =
                min(scalar(1), localFactor_[celli]*localRecovery_);
        }
    }

    // Base factors of this iteration's trials; no limit event yet
    localFactor0_ = localFactor_;
    localLimitedNow_ = false;
}


Foam::label Foam::ptcControl::applyLocalLimit
(
    scalarField& rDeltaTV,
    const vectorField& rMom,
    const scalarField& aMom,
    const scalar Uref,
    scalarField* ratio
)
{
    nLocalLimited_ = 0;
    nLocalThrottled_ = 0;
    nLocalSticky_ = 0;

    if (!localLimit_)
    {
        if (ratio)
        {
            *ratio = Zero;
        }
        return 0;
    }

    const scalar dUmax = fLoc_*Uref;

    if (!localMemory_)
    {
        // Memoryless limiter (7.3); implicit: |r| / (V/dt + a_P)
        forAll(rDeltaTV, celli)
        {
            // dU = |r| dt/V = |r| / (V/dt); GUARD: rDeltaTV >= cfVSmall/CFL
            const scalar dU =
                mag(rMom[celli])
               /max
                (
                    rDeltaTV[celli] + (localImplicit_ ? aMom[celli] : 0),
                    cfVSmall<scalar>()
                );
            if (ratio)
            {
                (*ratio)[celli] = dU/max(dUmax, cfVSmall<scalar>());
            }
            if (dU > dUmax)
            {
                if (localImplicit_)
                {
                    rDeltaTV[celli] =
                        mag(rMom[celli])/max(dUmax, cfVSmall<scalar>()) - aMom[celli];
                }
                else
                {
                    // dt <- dt*dUmax/dU  <=>  V/dt <- V/dt * dU/dUmax
                    rDeltaTV[celli] *= dU/max(dUmax, cfVSmall<scalar>());
                }
                ++nLocalLimited_;
            }
        }

        reduce(nLocalLimited_, sumOp<label>());
        return nLocalLimited_;
    }

    // Limiter with memory (D-055): dt_P <- f_P dt_P first, then the check.
    // f_P is the factor at the start of the iteration, so a line-search
    // retrial (smaller global CFL) does not compound the cut of an earlier
    // trial; the last trial's factor is the one kept. A cell counts one
    // limit event per outer iteration, however many trials limit it.
    forAll(rDeltaTV, celli)
    {
        // GUARD: f_P in (0, 1] by construction
        const scalar f = max(localFactor0_[celli], cfVSmall<scalar>());
        localFactor_[celli] = f;
        rDeltaTV[celli] /= f;

        const scalar aP = (localImplicit_ ? aMom[celli] : 0);
        const scalar dU = mag(rMom[celli])/max(rDeltaTV[celli] + aP, cfVSmall<scalar>());
        if (ratio)
        {
            (*ratio)[celli] = dU/max(dUmax, cfVSmall<scalar>());
        }
        if (dU > dUmax)
        {
            // Explicit: V/dt <- V/dt dU/dUmax. Implicit: V/dt <- the value
            // that makes |r|/(V/dt + a_P) = dUmax, i.e. |r|/dUmax - a_P
            // (> V/dt since dU > dUmax)
            const scalar rDT0 = rDeltaTV[celli];
            rDeltaTV[celli] = mag(rMom[celli])/max(dUmax, cfVSmall<scalar>()) - aP;
            const scalar cut = rDT0/max(rDeltaTV[celli], cfVSmall<scalar>());
            localFactor_[celli] = f*cut;
            localHoldRemaining_[celli] = localHold_;
            if (!localLimitedNow_[celli])
            {
                localLimitedNow_[celli] = true;
                ++localCount_[celli];
            }
            ++nLocalLimited_;
        }
        if (localFactor_[celli] < 1)
        {
            ++nLocalThrottled_;
        }
        if (localStickyAfter_ > 0 && localCount_[celli] >= localStickyAfter_)
        {
            ++nLocalSticky_;
        }
    }

    reduce(nLocalLimited_, sumOp<label>());
    reduce(nLocalThrottled_, sumOp<label>());
    reduce(nLocalSticky_, sumOp<label>());
    return nLocalLimited_;
}


void Foam::ptcControl::update(const doubleScalar R)
{
    const doubleScalar Rold = Rprev_;
    Rprev_ = R;
    // The CFL this accepted iteration was taken at; recordCost() pairs it
    // with the iteration's wall time (EFF, D-083)
    pendingCFL_ = CFL_;

    if (Rold <= 0)
    {
        return;
    }

    // GUARD: R > 0 before division
    const doubleScalar ratio = Rold/max(R, cfVSmall<doubleScalar>());

    doubleScalar CFLnew = CFL_;
    word why;

    switch (strategy_)
    {
        case strategy::mRDM:
        {
            if (R <= Rold)
            {
                // GUARD: pow base > 0 (ratio >= 1 here)
                const doubleScalar f = min(betaMax_, std::pow(ratio, gamma_));
                CFLnew = min(CFLmax_, CFL_*max(doubleScalar(1), f));
                why = "up:residualFell";
            }
            else
            {
                why = "noUp:residualRose";
            }
            break;
        }
        case strategy::EXP:
        {
            CFLnew = min(CFLmax_, CFL_*betaExp_);
            why = "up:EXP";
            break;
        }
        case strategy::SER:
        {
            CFLnew = min(CFLmax_, max(CFLmin_, CFL_*ratio));
            why = (ratio >= 1 ? "up:SER" : "down:SERresidualRose");
            break;
        }
        case strategy::EFF:
        {
            // D-083: the CFL of the best turnaround, from the share of the
            // linear solve in the wall time of an iteration (a ratio within
            // one iteration: most of the scatter of the wall time cancels).
            // If the solve cost grows like CFL^b, the pseudo-time per second
            // is largest at a share of 1/b; on the F1 case the optimum (CFL
            // ~43, 3.7 CFL/s) sat at 62 %. Hill climbing on CFL/seconds,
            // the first version, failed on the +-25 % scatter of the wall
            // time at low CFL.
            //
            // Decide on a window, not on one iteration (user): after a CFL
            // change skip one iteration, then average effWindow iterations
            // at that CFL. The step aims at the target share,
            //   ts ~ CFL^b,  s = ts/(t0 + ts)
            //   => f = [(s*/(1 - s*))/(s/(1 - s))]^(1/b),
            // clamped to 0.5..2 per decision: large steps far from the
            // target, small ones near it.
            if (!haveSample_)
            {
                why = "keep:noCleanSample";
                break;
            }
            haveSample_ = false;
            // GUARD: sampleCost_ > 0 (a measured wall time)
            const doubleScalar c = max(sampleCost_, cfVSmall<doubleScalar>());
            effLast_ = sampleCFL_/c;

            if (winCFL_ <= 0 || mag(sampleCFL_ - winCFL_) > 1e-6*sampleCFL_)
            {
                // First sample at a new CFL: the transient of the change
                winCFL_ = sampleCFL_;
                winN_ = 0;
                winShare_ = 0;
                winCost_ = 0;
                why = "keep:measuring";
                break;
            }
            winShare_ += sampleSolve_/c;
            winCost_ += c;
            ++winN_;
            if (winN_ < effWindow_)
            {
                why = "keep:measuring";
                break;
            }
            // Window complete
            const doubleScalar s =
                min(max(winShare_/winN_, doubleScalar(1e-3)), doubleScalar(0.999));
            shareSm_ = s;
            effLast_ = winCFL_/(winCost_/winN_);
            winCFL_ = -1;               // the next sample starts a window

            const doubleScalar lo = effShare_*(1 - effTol_);
            const doubleScalar hi = effShare_*(1 + effTol_);
            if (s >= lo && s <= hi)
            {
                why = "keep:solverShareOnTarget";
                break;
            }
            const doubleScalar f = min
            (
                max
                (
                    std::pow
                    (
                        (effShare_/(1 - effShare_))/(s/(1 - s)),
                        1/effExponent_
                    ),
                    doubleScalar(0.5)
                ),
                doubleScalar(2)
            );
            if (f > 1)
            {
                if (R > 2*Rold)
                {
                    why = "noUp:residualJump";
                    break;
                }
                CFLnew = min(CFLmax_, CFL_*f);
                why = "up:solverShareBelowTarget";
            }
            else
            {
                // The linear solve costs more than the optimum allows: the
                // reason is solver cost, not an instability
                CFLnew = max(CFLmin_, CFL_*f);
                why = "down:solverShareAboveTarget(solverCost)";
            }
            break;
        }
    }

    if (holdRemaining_ > 0)
    {
        --holdRemaining_;
        // No increase during the hold (SER may still decrease)
        if (CFLnew > CFL_)
        {
            why = "hold:afterCut";
        }
        CFLnew = min(CFLnew, CFL_);
    }

    // D-081: the ceiling from the last failed linear solve, relaxed by one
    // step per accepted iteration so that the CFL can probe upwards again
    if (failCeiling_ > 0)
    {
        if (CFLnew > max(ceiling_, CFL_))
        {
            why = "capped:failCeiling";
        }
        CFLnew = min(CFLnew, max(ceiling_, CFL_));
        ceiling_ = min(CFLmax_, ceiling_*ceilingRelax_);
    }
    if (CFLnew >= CFLmax_ && CFL_ < CFLmax_)
    {
        why = "capped:CFLmax";
    }

    if (!why.empty())
    {
        reason_ = (reason_.empty() ? why : word(reason_ + "+" + why));
    }
    CFL_ = CFLnew;
}


void Foam::ptcControl::recordCost
(
    const doubleScalar seconds,
    const label linIters,
    const doubleScalar solveSeconds
)
{
    if (pendingCFL_ > 0 && seconds > 0)
    {
        sampleCFL_ = pendingCFL_;
        sampleCost_ = seconds;
        sampleLin_ = linIters;
        sampleSolve_ = solveSeconds;
        haveSample_ = true;
    }
}


void Foam::ptcControl::linearFailure()
{
    if (failCeiling_ > 0)
    {
        ceiling_ = max(CFLmin_, failCeiling_*CFL_);
    }
}


void Foam::ptcControl::decrease(const doubleScalar factor)
{
    CFL_ = max(CFLmin_, factor*CFL_);
    holdRemaining_ = nHold_;
    // EFF (D-083): a cut is a restart of the search, not a measurement.
    // Without this the next sample (lower CFL, same cost) compared worse
    // than the one before the cut and was logged as solver cost - on the
    // F1 half-car right after a sentinel rollback, with 1 linear iteration.
    effPrev_ = -1;
    effPrevCFL_ = -1;
    effBest_ = -1;
    shareSm_ = -1;
    winCFL_ = -1;
    effDir_ = 1;
    effStreak_ = 0;
    haveSample_ = false;
}


void Foam::ptcControl::decrease(const doubleScalar factor, const word& why)
{
    decrease(factor);
    reason_ = (reason_.empty() ? why : word(reason_ + "+" + why));
}


void Foam::ptcControl::boost(const doubleScalar factor)
{
    if (holdRemaining_ > 0 || !(factor > 1))
    {
        return;
    }
    CFL_ = min(CFLmax_, factor*CFL_);
    if (failCeiling_ > 0)
    {
        CFL_ = min(CFL_, max(ceiling_, CFLmin_));     // D-081
    }
}


void Foam::ptcControl::writeState(dictionary& dict) const
{
    dict.set("CFL", CFL_);
    dict.set("Rprev", Rprev_);
    dict.set("nHoldRemaining", holdRemaining_);
    dict.set("cflCeiling", ceiling_);
    dict.set("effPrev", effPrev_);
    dict.set("effBest", effBest_);
    dict.set("effDir", effDir_);
    dict.set("effStep", effStep_);
    dict.set("effStreak", effStreak_);

    if (localMemory_)
    {
        // Sparse: cells with any memory (local labels, D-034)
        DynamicList<label> cells, holds, counts;
        DynamicList<scalar> factors;
        forAll(localFactor_, celli)
        {
            if
            (
                localFactor_[celli] < 1
             || localHoldRemaining_[celli] > 0
             || localCount_[celli] > 0
            )
            {
                cells.append(celli);
                factors.append(localFactor_[celli]);
                holds.append(localHoldRemaining_[celli]);
                counts.append(localCount_[celli]);
            }
        }
        dict.set("localLimitCells", labelList(cells));
        dict.set("localLimitFactor", scalarList(factors));
        dict.set("localLimitHold", labelList(holds));
        dict.set("localLimitCount", labelList(counts));
    }
}


void Foam::ptcControl::readState(const dictionary& dict)
{
    CFL_ = dict.get<doubleScalar>("CFL");
    Rprev_ = dict.getOrDefault<doubleScalar>("Rprev", -1);
    holdRemaining_ = dict.getOrDefault<label>("nHoldRemaining", 0);
    ceiling_ = dict.getOrDefault<doubleScalar>("cflCeiling", CFLmax_);
    effPrev_ = dict.getOrDefault<doubleScalar>("effPrev", -1);
    effBest_ = dict.getOrDefault<doubleScalar>("effBest", -1);
    effDir_ = dict.getOrDefault<label>("effDir", 1);
    effStep_ = dict.getOrDefault<doubleScalar>("effStep", betaExp_);
    effStreak_ = dict.getOrDefault<label>("effStreak", 0);

    if (localMemory_)
    {
        localFactor_ = scalar(1);
        localHoldRemaining_ = 0;
        localCount_ = 0;

        const labelList cells
        (
            dict.getOrDefault<labelList>("localLimitCells", labelList())
        );
        const scalarList factors
        (
            dict.getOrDefault<scalarList>("localLimitFactor", scalarList())
        );
        const labelList holds
        (
            dict.getOrDefault<labelList>("localLimitHold", labelList())
        );
        const labelList counts
        (
            dict.getOrDefault<labelList>("localLimitCount", labelList())
        );

        if
        (
            cells.size() != factors.size()
         || cells.size() != holds.size()
         || cells.size() != counts.size()
        )
        {
            FatalIOErrorInFunction(dict)
                << "localLimitCells/Factor/Hold/Count differ in size"
                << exit(FatalIOError);
        }

        forAll(cells, i)
        {
            if (cells[i] < 0 || cells[i] >= mesh_.nCells())
            {
                FatalIOErrorInFunction(dict)
                    << "localLimitCells cell " << cells[i]
                    << " out of range (restart with a different"
                    << " decomposition?)" << exit(FatalIOError);
            }
            if (!(factors[i] > 0) || factors[i] > 1)
            {
                FatalIOErrorInFunction(dict)
                    << "localLimitFactor " << factors[i]
                    << " outside (0, 1]" << exit(FatalIOError);
            }
            localFactor_[cells[i]] = factors[i];
            localHoldRemaining_[cells[i]] = holds[i];
            localCount_[cells[i]] = counts[i];
        }
    }
}


void Foam::ptcControl::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("cflStrategy", strategyName_);
    d.add("CFL0", CFL0_);
    d.add("CFLmin", CFLmin_);
    d.add("CFLmax", CFLmax_);
    d.add("gamma", gamma_);
    d.add("betaMax", betaMax_);
    d.add("betaExp", betaExp_);
    d.add("nHold", nHold_);
    d.add("continuityFactor", continuityFactor_);
    d.add("failCeiling", failCeiling_);
    d.add("ceilingRelax", ceilingRelax_);
    d.add("effTol", effTol_);
    d.add("effSolveShare", effShare_);
    d.add("effWindow", effWindow_);
    d.add("effExponent", effExponent_);
    dict.add("ptc", d);

    dictionary l;
    l.add("enabled", localLimit_);
    l.add("fLoc", fLoc_);
    l.add("implicit", localImplicit_);
    l.add("memory", localMemory_);
    l.add("localRecovery", localRecovery_);
    l.add("localHold", localHold_);
    l.add("localStickyAfter", localStickyAfter_);
    dict.add("localLimit", l);
}


// ************************************************************************* //
