/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    coupledFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state, incompressible, turbulent flow solver with block-coupled
    pressure-velocity solution (4x4 blocks u, v, w, p per cell), single-
    precision block linear algebra with double-precision reductions, block-
    GAMG-preconditioned Krylov solvers, pseudo-transient continuation with
    adaptive CFL, physicality line search, remediation cell sets, implicit
    MRF, FPE safety with sentinel and rollback, and exact restart.
    Turbulence is solved segregated after every coupled update.

    Specification: SPEC_coupledFoam.md (v2). Deviations: DECISIONS.md.

    Outer iteration (spec 5.7):
        nuEff, local dt (5.4) with the current CFL and remediation factors
        assemble momentum rows, local CFL limit (7.3), assemble continuity
        solve A dx = b - A x
        line search omega (7.2); if omega < omegaMin: cut CFL, repeat
        U += omega dU, p += omega dp, flux update (5.3e)
        sentinel (9.3), dynamic remediation (8.2)
        turbulence->correct() (5.8), bounds
        CFL strategy update (7.1), start-up switch (7.4), log (12.1)
        write (fields, coupledState), convergence (12.3)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "bound.H"
#include "clockTime.H"

#include "coupledDefaults.H"
#include "blockScalar.H"
#include "blockSolver.H"
#include "doubleReduce.H"
#include "blockGAMGPrecon.H"
#include "coupledAssembler.H"
#include "MRFCoupling.H"
#include "ptcControl.H"
#include "lineSearch.H"
#include "remediation.H"
#include "sentinel.H"
#include "convergenceMonitor.H"
#include "coupledState.H"
#include "runInfo.H"
#include "jsonWriter.H"
#include "gamgAutoTune.H"
#include "anderson.H"
#include "adaptiveTolerance.H"
#include "startupControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady incompressible solver with block-coupled pressure-velocity"
        " solution (coupledFoam)."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * Controls * * * * * * * * * * * * * * * * //

    const dictionary coupledDict(mesh.solutionDict().subOrEmptyDict("coupled"));
    const dictionary& linearDict = mesh.solverDict("coupled");

    const label maxIter =
        coupledDict.getOrDefault<label>("maxIter", coupledDefaults::outerMaxIter);
    const bool potentialInit =
        coupledDict.getOrDefault<bool>
        (
            "potentialInit",
            coupledDefaults::potentialInit
        );
    const bool ftz =
        coupledDict.getOrDefault<bool>("ftz", coupledDefaults::ftz);

    const dictionary& boundsDict = coupledDict.subOrEmptyDict("bounds");
    const scalar kMin =
        boundsDict.getOrDefault<scalar>("kMin", coupledDefaults::kMin);
    const scalar omegaMinBound =
        boundsDict.getOrDefault<scalar>
        (
            "omegaMin",
            coupledDefaults::boundOmegaMin
        );
    const scalar nutMaxFactor =
        boundsDict.getOrDefault<scalar>
        (
            "nutMaxFactor",
            coupledDefaults::nutMaxFactor
        );

    // FTZ/DAZ for benchmark runs (spec 9.1)
    const bool ftzApplied = (ftz ? runInfo::enableFTZ() : false);

    // * * * * * * * * * * * * * Components * * * * * * * * * * * * * * * * //

    MRFCoupling mrfc(mesh, MRF);
    coupledAssembler assembler(mesh, coupledDict, p, &mrfc);
    autoPtr<blockSolver> linSolver =
        blockSolver::New(assembler.matrix(), linearDict);

    ptcControl ptc(mesh, coupledDict);
    lineSearch ls(coupledDict);
    remediation rem(mesh, coupledDict);
    sentinel sen(mesh, coupledDict);
    anderson aa(mesh, coupledDict);
    adaptiveTolerance ew(linearDict);
    startupControl startup(coupledDict);

    // Linear-system dump for offline preconditioner studies
    // (Test-blockSystem): coupled.dumpLinearSystem (iterations), serial only
    const labelList dumpIters
    (
        coupledDict.getOrDefault<labelList>("dumpLinearSystem", labelList())
    );
    const label maxLinFails =
        coupledDict.subOrEmptyDict("ptc").getOrDefault<label>
        (
            "maxLinFails",
            coupledDefaults::maxLinFails
        );
    // B4 failure definition (D-049): strict = not converged to eta within
    // maxIter; reduction = a capped solve is accepted if it is finite and
    // reduced the true residual to at most linAcceptReduction times the
    // initial one
    const word linFailPolicy =
        coupledDict.subOrEmptyDict("ptc").getOrDefault<word>
        (
            "linFailPolicy",
            word(coupledDefaults::linFailPolicy)
        );
    const scalar linAcceptReduction =
        coupledDict.subOrEmptyDict("ptc").getOrDefault<scalar>
        (
            "linAcceptReduction",
            coupledDefaults::linAcceptReduction
        );
    if (linFailPolicy != "strict" && linFailPolicy != "reduction")
    {
        FatalIOErrorInFunction(coupledDict)
            << "ptc.linFailPolicy must be strict or reduction, got "
            << linFailPolicy << exit(FatalIOError);
    }
    if (!(linAcceptReduction > 0) || !(linAcceptReduction < 1))
    {
        FatalIOErrorInFunction(coupledDict)
            << "ptc.linAcceptReduction must be in (0, 1), got "
            << linAcceptReduction << exit(FatalIOError);
    }
    const bool linFailReduction = (linFailPolicy == "reduction");
    label nLinAccepted = 0;

    if (maxLinFails < 1)
    {
        FatalIOErrorInFunction(coupledDict)
            << "ptc.maxLinFails must be >= 1, got " << maxLinFails
            << exit(FatalIOError);
    }

    // autoTune controller of the block-GAMG cycle (6.3.5)
    autoPtr<gamgAutoTune> tuner;
    {
        const blockGAMGPrecon* gp =
            dynamic_cast<const blockGAMGPrecon*>(linSolver->preconditioner());
        if (gp)
        {
            tuner.reset
            (
                new gamgAutoTune
                (
                    gp->gamg(),
                    linearDict.subOrEmptyDict("blockGAMG")
                )
            );
        }
    }
    convergenceMonitor conv(coupledDict);
    coupledState state(mesh, coupledDict);

    volScalarField* kPtr = mesh.getObjectPtr<volScalarField>("k");
    volScalarField* omegaPtr = mesh.getObjectPtr<volScalarField>("omega");
    volScalarField* nutPtr = mesh.getObjectPtr<volScalarField>("nut");

    // * * * * * * * * * * * Effective settings (11) * * * * * * * * * * * * //

    {
        dictionary eff;
        eff.add("maxIter", maxIter);
        eff.add("potentialInit", potentialInit);
        startup.writeSettings(eff);
        eff.add("nonOrthLimiter", assembler.noc().limiter());
        ptc.writeSettings(eff);
        ls.writeSettings(eff);
        rem.writeSettings(eff);
        eff.subDict("remediation").subDict("static").add
        (
            "nonOrthLimiter",
            assembler.noc().limiterStatic()
        );
        sen.writeSettings(eff);
        {
            dictionary b;
            b.add("kMin", kMin);
            b.add("omegaMin", omegaMinBound);
            b.add("nutMaxFactor", nutMaxFactor);
            eff.add("bounds", b);
        }
        conv.writeSettings(eff);
        eff.add("ftz", ftz);
        eff.add("writeState", state.writeState());
        eff.add
        (
            "guards",
            dictionary(coupledDict.subOrEmptyDict("guards"))
        );
        dictionary lin;
        linSolver->writeSettings(lin);
        ew.writeSettings(lin);
        eff.add("linearSolver", lin);
        {
            dictionary a;
            aa.writeSettings(a);
            eff.add("anderson", a);
        }
        eff.subDict("ptc").add("maxLinFails", maxLinFails);
        eff.subDict("ptc").add("linFailPolicy", linFailPolicy);
        eff.subDict("ptc").add("linAcceptReduction", linAcceptReduction);
        if (tuner)
        {
            eff.add("autoTune", tuner->settings());
        }

        Info<< nl << "coupledFoam: effective settings" << nl
            << eff << endl;
    }

    // * * * * * * * * * * * * Once-per-run block (12.2) * * * * * * * * * * //

    Info<< "coupledFoam: mesh cells "
        << returnReduce(mesh.nCells(), sumOp<label>())
        << ", internal faces "
        << returnReduce(mesh.nInternalFaces(), sumOp<label>())
        << ", ranks " << UPstream::nProcs() << nl
        << "coupledFoam: precision " << runInfo::precisionDict() << nl
        << "coupledFoam: FTZ " << (ftzApplied ? "on" : "off")
        << ", FPE traps " << (runInfo::fpeActive() ? "on" : "off")
        << ", MRF " << (mrfc.active() ? "active" : "none") << nl
        << "coupledFoam: patches:";
    forAll(assembler.kinds(), patchi)
    {
        Info<< ' ' << mesh.boundary()[patchi].name() << '='
            << boundaryCoupling::kindName(assembler.kinds()[patchi]);
    }
    Info<< endl;

    rem.buildStatic();
    assembler.setStaticCells(rem.isStatic());

    // * * * * * * * * * * * * * * * Restart (10) * * * * * * * * * * * * * //

    label iter = 0;
    scalar R1 = -1;
    scalar phiConsistency = -1;
    bool restarted = false;

    {
        dictionary st;
        if (state.read(st))
        {
            restarted = true;
            iter = st.get<label>("iter");
            startup.readState(st, iter);
            R1 = st.get<scalar>("R1");
            ls.setReference(st.get<scalar>("Uref"), st.get<scalar>("pref"));
            ptc.readState(st);
            rem.readState(st);
            sen.readState(st);
            ew.readState(st);
            conv.readState(st);
            if (tuner)
            {
                tuner->readState(st);
            }

            Info<< "coupledFoam: restart from " << runTime.timeName()
                << " at iteration " << iter << ", CFL " << ptc.CFL()
                << ", start-up beta " << startup.beta(iter + 1)
                << " (trigger "
                << (startup.trigger().empty() ? word("pending") : startup.trigger())
                << ")" << endl;

            // phi consistency with the stored D (spec 10)
            tmp<volScalarField> tD = state.readD();
            if (tD.valid())
            {
                const scalarField& Dc = tD().primitiveField();
                scalarField abar(mesh.V()/max(Dc, VSMALL));  // GUARD
                assembler.rc().updateD(abar, p);

                // The written phi was built with the explicit term q_f of
                // the assembly that was solved (rhieChow::updateFlux); use
                // the stored q_f if present, else recompute it from p
                surfaceScalarField phiRe("phiRecomputed", phi);
                tmp<surfaceScalarField> tQ = state.readQ();
                if (tQ.valid())
                {
                    assembler.rc().setQ(tQ());
                }
                else
                {
                    Info<< "coupledFoam: restart without coupledQ, q_f"
                        << " recomputed from p for the phi check" << endl;
                }
                assembler.rc().updateFlux
                (
                    phiRe, U, p, assembler.noc(), !tQ.valid()
                );
                MRF.makeRelative(phiRe);

                scalar dmax = gMax(mag(phi.primitiveField() - phiRe.primitiveField())());
                scalar pmax = gMax(mag(phi.primitiveField())());
                forAll(phi.boundaryField(), patchi)
                {
                    const scalarField& a = phi.boundaryField()[patchi];
                    const scalarField& b = phiRe.boundaryField()[patchi];
                    if (a.size())
                    {
                        dmax = max(dmax, max(mag(a - b)()));
                        pmax = max(pmax, max(mag(a)()));
                    }
                }
                reduce(dmax, maxOp<scalar>());
                reduce(pmax, maxOp<scalar>());
                // GUARD
                phiConsistency = dmax/max(pmax, VSMALL);
                Info<< "coupledFoam: restart phi consistency "
                    << phiConsistency << endl;
            }
            else
            {
                Info<< "coupledFoam: restart without coupledD, phi"
                    << " consistency not checked" << endl;
            }
        }
        else
        {
            Info<< "coupledFoam: fresh start"
                << (potentialInit ? " (potentialInit: fields from"
                    " potentialFoam expected)" : "") << endl;
        }
    }

    // * * * * * * * * * * * * * * * Helpers * * * * * * * * * * * * * * * * //

    const Vector<label>& solD = mesh.solutionD();

    auto stateDict = [&]() -> dictionary
    {
        dictionary st;
        st.set("iter", iter);
        st.set("startupDone", startup.done(iter + 1));
        startup.writeState(st);
        st.set("R1", R1);
        st.set("Uref", ls.Uref());
        st.set("pref", ls.pref());
        ptc.writeState(st);
        rem.writeState(st);
        sen.writeState(st);
        ew.writeState(st);
        conv.writeState(st);
        if (tuner)
        {
            tuner->writeState(st);
        }
        // Anderson history deliberately not part of the state (D-026)
        st.set("refinementHistory", labelList());
        return st;
    };

    auto writeOutputs = [&]()
    {
        rem.write();
        state.writeD(assembler.rc().D());
        state.writeQ(assembler.rc().q());
        state.write(stateDict());
    };

    auto applyBounds = [&]() -> label
    {
        if (kPtr)
        {
            bound(*kPtr, dimensionedScalar(kPtr->dimensions(), kMin));
        }
        if (omegaPtr)
        {
            bound
            (
                *omegaPtr,
                dimensionedScalar(omegaPtr->dimensions(), omegaMinBound)
            );
        }
        label nCapped = 0;
        if (nutPtr)
        {
            const tmp<volScalarField> tnu = laminarTransport.nu();
            const scalarField& nu = tnu().primitiveField();
            scalarField& nut = nutPtr->primitiveFieldRef();
            forAll(nut, celli)
            {
                const scalar nutMax = nutMaxFactor*nu[celli];
                if (nut[celli] > nutMax)
                {
                    nut[celli] = nutMax;
                    ++nCapped;
                }
            }
            reduce(nCapped, sumOp<label>());
            if (nCapped)
            {
                nutPtr->correctBoundaryConditions();
            }
        }
        return nCapped;
    };

    label nCflCutsTotal = 0;
    label linFails = 0;
    bool converged = false;
    clockTime runTimer;
    scalar lastR = -1;
    label lastLinIters = 0;

    // * * * * * * * * * * * * * * * Outer loop * * * * * * * * * * * * * * * //

    Info<< "\ncoupledFoam: starting outer iterations\n" << endl;

    while (runTime.loop())
    {
        ++iter;
        clockTime iterTimer;
        scalar tAsm = 0, tSolve = 0, tTurb = 0;

        mrfc.update();
        MRF.correctBoundaryVelocity(U);

        if (!ls.referenceSet())
        {
            ls.setReference(U);
            Info<< "coupledFoam: Uref " << ls.Uref() << ", pref "
                << ls.pref() << endl;
        }

        const volScalarField nuEff("nuEff", turbulence->nuEff());
        // Start-up (D-048): developed-start probe on iteration 1 of a
        // fresh start without potentialInit and with a non-uniform U
        if (iter == 1 && !restarted)
        {
            const vectorField& Ui = U.primitiveField();
            const vector U0 = (Ui.size() ? Ui[0] : vector::zero);
            scalar dev = 0;
            forAll(Ui, celli)
            {
                dev = max(dev, mag(Ui[celli] - U0));
            }
            // Uniform on every rank and equal across ranks
            scalar devG = dev;
            reduce(devG, maxOp<scalar>());
            vector Umin = (Ui.size() ? U0 : vector::uniform(GREAT));
            vector Umax = (Ui.size() ? U0 : vector::uniform(-GREAT));
            reduce(Umin, minOp<vector>());
            reduce(Umax, maxOp<vector>());
            const bool nonUniform = devG > 0 || mag(Umax - Umin) > 0;
            startup.startProbe(!potentialInit && nonUniform);
        }
        scalar betaGlobal = startup.beta(iter);
        const bool startupDone = startup.done(iter);
        if (iter > 1 && betaGlobal != startup.beta(iter - 1))
        {
            // The discretisation changes along the ramp: the Anderson
            // history refers to another operator
            aa.flush();
        }
        scalarField beta(rem.beta(betaGlobal));
        const scalarField cflF(rem.cflFactor());

        blockScalarList dx(blockDim*mesh.nCells(), Zero);
        blockSolverPerformance perf;
        scalar omega = 1;
        label cuts = 0;
        label nLocLim = 0;
        scalar Rraw = 0;
        scalar eta = 0;
        bool skipStep = false;

        // --- Assemble, solve, line search with CFL cuts (7.2)
        while (true)
        {
            clockTime ta;
            scalarField rDTV(ptc.rDeltaTV(phi, nuEff, cflF));
            assembler.assembleMomentum(U, p, phi, nuEff, beta);
            nLocLim = ptc.applyLocalLimit
            (
                rDTV,
                assembler.momentumResidual(),
                ls.Uref()
            );
            assembler.assembleContinuity(rDTV);
            tAsm += ta.elapsedTime();

            Rraw = assembler.residualL2();

            if (startup.probing())
            {
                const scalar rU1 = assembler.rU();
                const scalar rp1 = assembler.rp();
                if (startup.decideDeveloped(rU1, rp1))
                {
                    Info<< "coupledFoam: developed start (rU " << rU1
                        << ", rp " << rp1 << "): start-up skipped, beta 1"
                        << " from iteration 1" << endl;
                }
                else
                {
                    Info<< "coupledFoam: start not developed (rU " << rU1
                        << ", rp " << rp1 << "): iteration 1 re-assembled"
                        << " with the start-up beta" << endl;
                    betaGlobal = startup.beta(iter);
                    beta = rem.beta(betaGlobal);
                    continue;
                }
            }

            // Eisenstat-Walker inner tolerance (amendment B2)
            {
                // GUARD: R1 >= VSMALL (9.2)
                const scalar Rn = Rraw/max((R1 > 0 ? R1 : Rraw), VSMALL);
                eta = ew.eta(Rn, startupDone);
                linSolver->setRelTol(eta);
            }

            if (!UPstream::parRun() && dumpIters.found(iter) && cuts == 0)
            {
                const blockLduMatrix4& Am = assembler.matrix();
                const fileName dumpFile
                (
                    runTime.path()/"linsys"/("iter" + Foam::name(iter))
                );
                mkDir(dumpFile.path());
                OFstream os(dumpFile, IOstreamOption(IOstreamOption::BINARY));
                os  << Am.diag() << Am.upper() << Am.lower()
                    << assembler.rhs()
                    << token::SPACE << doubleScalar(assembler.normFactor())
                    << token::SPACE << doubleScalar(eta)
                    << token::SPACE << doubleScalar(ptc.CFL()) << nl;
                Info<< "coupledFoam: linear system written to " << dumpFile
                    << endl;
            }

            clockTime ts;
            dx = Zero;
            perf = linSolver->solve(dx, assembler.rhs(), assembler.normFactor());
            tSolve += ts.elapsedTime();

            // No increment in empty directions (2D)
            for (direction c = 0; c < vector::nComponents; ++c)
            {
                if (solD[c] == -1)
                {
                    for (label celli = 0; celli < mesh.nCells(); ++celli)
                    {
                        dx[celli*blockDim + c] = 0;
                    }
                }
            }

            // Linear-solve failure (amendment B4): not converged to
            // eta*||r0|| (or the absolute floor) within maxIter, or
            // non-finite values. CFL cut and repeat; maxLinFails
            // consecutive failures abort through the 9.3 path.
            const bool finiteDx = std::isfinite(doubleReduce::sumSqr(dx));
            const bool reduced =
                perf.finalResidual
             <= linAcceptReduction*perf.initialResidual;
            const bool solveFailed =
                !finiteDx
             || (linFailReduction ? !(perf.converged || reduced) : !perf.converged);
            if (!solveFailed && !perf.converged)
            {
                ++nLinAccepted;
            }

            if (solveFailed)
            {
                ++linFails;
                if (linFails >= maxLinFails)
                {
                    // Fields are unchanged since iteration iter-1: write
                    // them (and the remediation sets) into
                    // <iter-1>_lastValid, never into the current time
                    // directory, which would be left without U/p (9.3)
                    const label validIter = iter - 1;
                    const sentinel::checkResult chk =
                        sen.check(U, p, kPtr, omegaPtr, ls.Uref(), ls.pref());
                    sen.writeLastValid
                    (
                        U, p, phi, kPtr, omegaPtr, nutPtr, chk.offending,
                        validIter
                    );
                    rem.write(sentinel::lastValidName(validIter));
                    FatalErrorInFunction
                        << linFails << " consecutive linear-solve failures"
                        << " (maxLinFails " << maxLinFails << ") at iteration "
                        << iter << ", last: initial residual "
                        << perf.initialResidual << ", final "
                        << perf.finalResidual << " after "
                        << perf.nIterations << " iterations, eta " << eta
                        << ". Fields of iteration " << validIter
                        << " written to " << sentinel::lastValidName(validIter)
                        << "." << exit(FatalError);
                }
            }
            else
            {
                linFails = 0;
            }

            omega = (solveFailed ? 0 : ls.omega(dx));

            if
            (
                (solveFailed || omega < ls.omegaMin())
             && cuts < ls.maxCflCuts()
            )
            {
                ptc.decrease(ls.kappa());
                ++cuts;
                ++nCflCutsTotal;
                continue;
            }
            if (solveFailed)
            {
                // Cuts exhausted (maxCflCuts < maxLinFails): no update from
                // a failed solve, skip the step (D-020)
                skipStep = true;
                dx = Zero;
                omega = 0;
            }
            else if (omega < ls.omegaMin())
            {
                omega = ls.omegaMin();
                rem.markDynamic(ls.offendingCells(dx));
                aa.flush();
            }
            break;
        }

        // Anderson history is invalid after a CFL change (B5)
        if (cuts > 0 || skipStep)
        {
            aa.flush();
        }

        // GUARD: R1 >= VSMALL before division (9.2)
        if (R1 < 0)
        {
            R1 = max(Rraw, VSMALL);
        }
        const scalar R = Rraw/max(R1, VSMALL);
        const label nClamped = returnReduce(assembler.nClamped(), sumOp<label>());
        if (nClamped)
        {
            WarningInFunction
                << nClamped << " coefficients clamped at +-"
                << coupledDefaults::clampValue << " (spec 9.2)" << endl;
        }

        if (skipStep)
        {
            ptc.decrease(sen.cflFactor());
            Info<< "coupledFoam: linear solve failed after " << cuts
                << " CFL cuts at iteration " << iter
                << ", step skipped, CFL -> " << ptc.CFL() << endl;
        }

        // --- Field update
        if (!skipStep)
        {
            rem.clipIncrement(dx, U, omega, ls.Uref());
        }
        sen.store(U, p, phi, kPtr, omegaPtr, nutPtr);

        vectorField dUapplied(mesh.nCells());
        scalarField dpApplied(mesh.nCells());
        {
            vectorField& Ui = U.primitiveFieldRef();
            scalarField& pi = p.primitiveFieldRef();
            forAll(Ui, celli)
            {
                const blockScalar* d = dx.cdata() + celli*blockDim;
                dUapplied[celli] = omega*vector(d[0], d[1], d[2]);
                dpApplied[celli] = omega*scalar(d[blockP]);
                Ui[celli] += dUapplied[celli];
                pi[celli] += dpApplied[celli];
            }
        }
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // Anderson acceleration of the accepted update (amendment B5)
        label andersonStatus = -1;
        if (aa.enabled() && !skipStep)
        {
            const anderson::status st =
                aa.apply(U, p, dUapplied, dpApplied, ls.Uref(), ls.pref());
            andersonStatus = label(st);
            if (st == anderson::status::applied)
            {
                U.correctBoundaryConditions();
                p.correctBoundaryConditions();
                const sentinel::checkResult chkA =
                    sen.check(U, p, kPtr, omegaPtr, ls.Uref(), ls.pref());
                if (!chkA.ok)
                {
                    aa.reject(U, p);
                    U.correctBoundaryConditions();
                    p.correctBoundaryConditions();
                    andersonStatus = -2;
                }
            }
        }

        // Flux with the explicit Rhie-Chow term of the solved assembly: the
        // continuity row that was solved, evaluated at the new (U, p)
        // (conservative to the solver tolerance for omega = 1)
        assembler.rc().updateFlux(phi, U, p, assembler.noc(), false);
        MRF.makeRelative(phi);

        // --- Sentinel (9.3)
        bool rolledBack = false;
        auto rollback = [&](const sentinel::checkResult& chk)
        {
            sen.restore(U, p, phi, kPtr, omegaPtr, nutPtr);
            aa.flush();
            ptc.decrease(sen.cflFactor());
            rem.markDynamic(chk.offending);
            rolledBack = true;

            Info<< "coupledFoam: sentinel rollback at iteration " << iter
                << " (nonFinite " << chk.nNonFinite << ", max|U| "
                << chk.maxMagU << ", p [" << chk.minP << ", " << chk.maxP
                << "]), CFL -> " << ptc.CFL() << endl;

            if (sen.exhausted())
            {
                // The restored fields are those of iteration iter-1
                const label validIter = iter - 1;
                sen.writeLastValid
                (
                    U, p, phi, kPtr, omegaPtr, nutPtr, chk.offending,
                    validIter
                );
                rem.write(sentinel::lastValidName(validIter));
                FatalErrorInFunction
                    << sen.consecutive() << " consecutive rollbacks (limit "
                    << sen.maxRollbacks() << ") at iteration " << iter
                    << ". Last valid fields (iteration " << validIter
                    << ") written to " << sentinel::lastValidName(validIter)
                    << "." << exit(FatalError);
            }
        };

        {
            const sentinel::checkResult chk =
                sen.check(U, p, kPtr, omegaPtr, ls.Uref(), ls.pref());
            if (!chk.ok)
            {
                rollback(chk);
            }
        }

        label nNutCapped = 0;

        if (!rolledBack)
        {
            const label dynVersion = rem.dynamicVersion();
            rem.updateDynamic
            (
                U, p, kPtr, omegaPtr, ls.Uref(), ls.pref(), iter
            );
            if (rem.dynamicVersion() != dynVersion)
            {
                // Membership changed (beta and dt of those cells change)
                aa.flush();
            }

            // --- Turbulence, segregated (5.8)
            clockTime tt;
            laminarTransport.correct();
            turbulence->correct();
            nNutCapped = applyBounds();
            tTurb = tt.elapsedTime();

            const sentinel::checkResult chk =
                sen.check(U, p, kPtr, omegaPtr, ls.Uref(), ls.pref());
            if (!chk.ok)
            {
                rollback(chk);
            }
        }

        if (!rolledBack)
        {
            // A skipped step (failed solve, cuts exhausted) is neither an
            // accepted step nor a rollback: no strategy/EW/sentinel update
            if (!skipStep)
            {
                sen.accepted();
                ptc.update(R);
                ew.accept(R, eta);

                // Line-search beta (7.2): CFL boost on a full step
                if (omega >= 1 && cuts == 0)
                {
                    ptc.boost(ls.beta());
                }
            }

            if (startup.update(iter, R))
            {
                // beta 0 -> 1 changes the discretisation: the Anderson
                // history refers to the upwind operator
                aa.flush();
                Info<< "coupledFoam: start-up ramp starts at iteration "
                    << iter << " (trigger " << startup.trigger() << ", R "
                    << R << "), beta 1 from iteration "
                    << startup.rampEndIter() << endl;
            }
        }

        conv.record(runTime);
        lastR = R;
        lastLinIters = perf.nIterations;

        if (tuner)
        {
            // Only accepted, successful solves enter the rho window
            // (rho < 0: not recorded, the iteration still counts)
            const bool useRho = !rolledBack && !skipStep && perf.converged;
            if (tuner->record(iter, useRho ? scalar(perf.rho) : scalar(-1)))
            {
                // Different preconditioner from the next solve on
                aa.flush();
            }
        }

        // --- Log line (12.1)
        Info<< "CF| iter=" << iter
            << " CFL=" << ptc.CFL()
            << " omega=" << omega
            << " cuts=" << cuts
            << " R=" << R
            << " rU=" << assembler.rU()
            << " rp=" << assembler.rp()
            << " linIters=" << perf.nIterations
            << " linRes=" << perf.finalResidual
            << " tAsm=" << tAsm
            << " tSolve=" << tSolve
            << " tTurb=" << tTurb
            << " tIter=" << iterTimer.elapsedTime()
            << " tWall=" << runTime.elapsedClockTime()
            << " nStat=" << rem.nStatic()
            << " nDyn=" << rem.nDynamic()
            << " nLocLim=" << nLocLim
            << " nRollback=" << sen.nRollbacks()
            << " nClamped=" << nClamped
            << " eta=" << eta
            << " rho=" << perf.rho;
        if (aa.enabled())
        {
            Info<< " aa=" << andersonStatus;
        }
        if (nNutCapped)
        {
            Info<< " nNutCapped=" << nNutCapped;
        }
        if (conv.haveForces())
        {
            Info<< " Cd=" << conv.Cd() << " Cl=" << conv.Cl();
        }
        Info<< endl;

        // --- Write
        if (!rolledBack && conv.converged(R))
        {
            converged = true;
            Info<< "coupledFoam: converged at iteration " << iter
                << " (R " << R << ")" << endl;
            // As native solvers: write, end, and let runTime.loop() run the
            // function objects for the final state
            runTime.writeAndEnd();
            writeOutputs();
            continue;
        }

        if (iter >= maxIter)
        {
            Info<< "coupledFoam: maxIter " << maxIter << " reached" << endl;
            runTime.writeAndEnd();
            writeOutputs();
            continue;
        }

        if (runTime.writeTime())
        {
            runTime.write();
            writeOutputs();
        }
    }

    // * * * * * * * * * * * * * * Run summary * * * * * * * * * * * * * * * //

    const scalar cpuSum = returnReduce(scalar(runInfo::cpuSeconds()), sumOp<scalar>());
    const label rssKB = runInfo::peakRSSkB();
    const label rssMax = returnReduce(rssKB, maxOp<label>());
    const scalar rssSum = returnReduce(scalar(rssKB), sumOp<scalar>());

    Info<< nl << "coupledFoam: iterations " << iter
        << ", converged " << converged
        << ", final R " << lastR
        << ", CFL cuts " << nCflCutsTotal
        << ", rollbacks " << sen.nRollbacks() << nl
        << "coupledFoam: wall time " << runTimer.elapsedTime() << " s, CPU "
        << cpuSum << " s (" << cpuSum/3600.0 << " CPU-h, all ranks)" << nl
        << "coupledFoam: peak RSS max rank " << rssMax/1024.0 << " MB, sum "
        << rssSum/1024.0 << " MB" << endl;

    {
        jsonWriter j;
        j.add("solver", "coupledFoam");
        j.add("nProcs", UPstream::nProcs());
        j.add("nCells", returnReduce(mesh.nCells(), sumOp<label>()));
        j.add("iterations", iter);
        j.add("converged", converged);
        j.add("finalR", lastR);
        j.add("finalCFL", ptc.CFL());
        j.add("lastLinearIterations", lastLinIters);
        j.add("cflCuts", nCflCutsTotal);
        j.add("linFailPolicy", linFailPolicy);
        j.add("linAcceptedUnconverged", nLinAccepted);
        j.add("rollbacks", sen.nRollbacks());
        j.add("staticCells", rem.nStatic());
        j.add("dynamicCells", rem.nDynamic());
        j.add("wallSeconds", runTimer.elapsedTime());
        j.add("cpuSeconds", cpuSum);
        j.add("cpuHours", cpuSum/3600.0);
        j.add("peakRSS_MB_maxRank", rssMax/1024.0);
        j.add("peakRSS_MB_sum", rssSum/1024.0);
        j.add("restarted", restarted);
        j.add("startupMode", startup.modeName());
        j.add
        (
            "startupTrigger",
            startup.trigger().empty() ? word("pending") : startup.trigger()
        );
        j.add("rampStartIter", startup.rampStartIter());
        j.add("rampEndIter", startup.rampEndIter());
        j.add("phiConsistency", phiConsistency);
        j.add("ftz", ftzApplied);
        j.add("fpeTraps", runInfo::fpeActive());
        if (conv.haveForces())
        {
            j.add("Cd", conv.Cd());
            j.add("Cl", conv.Cl());
        }
        if (tuner)
        {
            j.addRaw("gamgTuneEvents", tuner->eventsJson());
            j.add("gamgTuneFailed", tuner->failed());
        }
        if (aa.enabled())
        {
            j.add("andersonApplied", aa.nApplied());
            j.add("andersonSkipped", aa.nSkipped());
            j.add("andersonRejected", aa.nRejected());
            j.add("andersonFlushed", aa.nFlushed());
        }
        const blockGAMGPrecon* gp =
            dynamic_cast<const blockGAMGPrecon*>(linSolver->preconditioner());
        if (gp)
        {
            const dictionary gs(gp->gamg().statsDict());
            j.add("gamgLevels", gs.get<label>("nLevels"));
            j.add("gamgMergeLevels", gs.get<label>("mergeLevels"));
            j.add("gamgCop", gs.get<scalar>("Cop"));
            j.addList("gamgCellsPerLevel", gs.get<labelList>("cellsPerLevel"));
        }
        j.write(runTime.globalPath()/"postProcessing"/"coupledFoam"/"summary.json");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
