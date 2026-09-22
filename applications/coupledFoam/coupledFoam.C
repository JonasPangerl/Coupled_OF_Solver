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
#include "diagnosticFields.H"
#include "sentinel.H"
#include "convergenceMonitor.H"
#include "coupledState.H"
#include "runInfo.H"
#include "jsonWriter.H"
#include "gamgAutoTune.H"
#include "anderson.H"
#include "sfdControl.H"
#include "adaptiveTolerance.H"
#include "startupControl.H"
#include "diagnostics.H"
#include "mixedFvPatchFields.H"
#include "SolverPerformance.H"
#include "Pair.H"
#include <algorithm>
#include <limits>
#include <vector>

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
    sfdControl sfd(mesh, coupledDict, ptc.nHold());

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
    // B4 failure definition (linFailPolicy): strict = not converged to eta within
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

    // Deep diagnostics (TASK 5, D-045): level 0 = off. Every hook below is
    // guarded with diag.active(n); at level 0 the solver and the assembler
    // get no diagnostics object and take no clock readings.
    diagnostics diag(coupledDict);
    // Write-time diagnostic fields (amendment C6)
    diagnosticFields dfields(mesh, coupledDict);
    if (diag.active(1))
    {
        linSolver->setDiagnostics(&diag);
        assembler.setTiming(true);
    }

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
        sfd.writeSettings(eff);
        eff.subDict("ptc").add("maxLinFails", maxLinFails);
        eff.subDict("ptc").add("linFailPolicy", linFailPolicy);
        eff.subDict("ptc").add("linAcceptReduction", linAcceptReduction);
        if (tuner)
        {
            eff.add("autoTune", tuner->settings());
        }
        {
            dictionary dg;
            dg.add("level", diag.level());
            dg.add("echo", diag.echo());
            eff.add("diagnostics", dg);
        }
        dfields.writeSettings(eff);

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
            sfd.readState(st);
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
                if (assembler.rc().tensorial())
                {
                    // Tensorial D_f (C2) from the stored coupledDT
                    tmp<volTensorField> tDT = state.readDT();
                    if (tDT.valid())
                    {
                        assembler.rc().setDT(tDT(), p);
                    }
                    else
                    {
                        Info<< "coupledFoam: restart without coupledDT,"
                            << " phi check uses the scalar D_f" << endl;
                    }
                }

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

    // * * * * * * * * * * * * * Diagnostics set-up * * * * * * * * * * * * * //

    const blockGAMGPrecon* gpDiag =
        dynamic_cast<const blockGAMGPrecon*>(linSolver->preconditioner());
    diagPhase phase;
    label diagHierVersion = -1;
    doubleScalar diagWritePrev = 0;

    // Level 3: inflow/outflow state of the faces of mixed-type patches
    // (inletOutlet, freestream*, ...), rank-local
    labelList bcPatches;
    List<boolList> bcOut;
    bool bcInit = false;

    if (diag.active(1))
    {
        diagJson h;
        h.beginObject();
        h.add("type", "header");
        h.add("format", label(1));
        h.add("solver", "coupledFoam");
        h.add("level", diag.level());
        h.add("rank", label(UPstream::myProcNo()));
        h.add("nProcs", label(UPstream::nProcs()));
        h.add("nCellsLocal", label(mesh.nCells()));
        h.add("restarted", restarted);
        h.add("startIter", iter);
        h.add("case", std::string(runTime.globalPath()));
        // Fields evaluated on this rank only (no reduction added, TASK
        // 5.4): merge over the rank files (bench/diag_tools.py)
        h.beginArray("localKeys");
        for
        (
            const char* k
          : {
                "residuals.massErrMax", "residuals.massErrSum",
                "controls.dt", "turbulence.nBoundK",
                "turbulence.nBoundOmega", "timings", "memory",
                "controls.trials.violU", "controls.trials.violP",
                "bcFlips", "gamgSetup", "localLimit", "dynamicSet",
                "controls.sfd.maxDev"
            }
        )
        {
            h.valueRaw("\"" + std::string(k) + "\"");
        }
        h.endArray();
        h.endObject();
        diag.open(runTime, restarted, h.str());

        if (diag.active(3))
        {
            DynamicList<label> pl;
            forAll(mesh.boundary(), patchi)
            {
                if (mesh.boundary()[patchi].coupled())
                {
                    continue;
                }
                if
                (
                    isA<mixedFvPatchVectorField>(U.boundaryField()[patchi])
                 || isA<mixedFvPatchScalarField>(p.boundaryField()[patchi])
                )
                {
                    pl.append(patchi);
                }
            }
            bcPatches.transfer(pl);
            bcOut.resize(bcPatches.size());
        }

        Info<< "coupledFoam: diagnostics level " << diag.level()
            << " -> " << runTime.globalPath()/"diagnostics" << endl;
    }

    // * * * * * * * * * * * * * * * Helpers * * * * * * * * * * * * * * * * //

    const Vector<label>& solD = mesh.solutionD();

    auto stateDict = [&]() -> dictionary
    {
        // Entries are formatted with IOstream::defaultPrecision() when they
        // are created (6 digits by default): round-trip precision here, or
        // the restart scalars (CFL, R1, Uref, ...) are quantised before
        // coupledState writes them with max_digits10 (D-034, D-045)
        const unsigned oldPrecision = IOstream::defaultPrecision
        (
            std::numeric_limits<doubleScalar>::max_digits10
        );
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
        sfd.writeState(st);
        conv.writeState(st);
        if (tuner)
        {
            tuner->writeState(st);
        }
        // Anderson history deliberately not part of the state (D-026)
        st.set("refinementHistory", labelList());
        IOstream::defaultPrecision(oldPrecision);
        return st;
    };

    auto writeOutputs = [&]()
    {
        rem.write();
        // C6 (C3 hook: dfields.setSFD(&Ubar) while SFD is active -> USFD)
        dfields.write(phi);
        state.writeD(assembler.rc().D());
        if (assembler.rc().tensorial())
        {
            state.writeDT(assembler.rc().DT());
        }
        state.writeQ(assembler.rc().q());
        sfd.write();
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

    // Anderson history flush (7.5); the SFD filter is reset with it (C3)
    auto flushHistory = [&](const char* reason)
    {
        aa.flush(reason);
        sfd.reset(U);
    };

    // Force-coefficient window statistics per iteration (summary.json)
    DynamicList<label> forceHistIter;
    List<DynamicList<doubleScalar>> forceHistMean(3);
    List<DynamicList<doubleScalar>> forceHistRms(3);
    List<DynamicList<doubleScalar>> forceHistDrift(3);

    label nCflCutsTotal = 0;
    label nPivotFallbackTotal = 0;
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

        // Diagnostics state of this iteration (filled only if active)
        const scalar CFLstart = ptc.CFL();
        const label nFlushedStart = aa.nFlushed();
        label nSenChecks = 0;
        doubleScalar tFlux = 0;
        doubleScalar tDiagMain = 0;
        doubleScalar dtMin = 0, dtMed = 0, dtMax = 0;
        doubleScalar massErrMax = 0, massErrSum = 0;
        label nBoundK = 0, nBoundOmega = 0;
        bool turbSolved = false;
        diagJson trialsJ;
        if (diag.active(1))
        {
            diag.beginIteration();
            assembler.resetTimes();
            trialsJ.beginArray();
        }

        mrfc.update();
        MRF.correctBoundaryVelocity(U);

        if (!ls.referenceSet())
        {
            ls.setReference(U);
            Info<< "coupledFoam: Uref " << ls.Uref() << " (" << ls.UrefSource()
                << ", mode " << ls.UrefMode() << "), pref " << ls.pref()
                << endl;
        }
        if (sfd.enabled() && sfd.DeltaStar() == 0)
        {
            sfd.setReference(ls.Uref());
        }
        // SFD (7.6): activation after the start-up phase; C6 writes USFD
        // while it is active
        sfd.begin(U, startup.done(iter));
        assembler.setSFD
        (
            sfd.chiStar(),
            (sfd.active() ? &sfd.Ubar().primitiveField() : nullptr)
        );
        dfields.setSFD(sfd.active() ? &sfd.Ubar() : nullptr);

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
            reduce(dev, maxOp<scalar>());
            // Uniform on every rank, but different values across ranks
            vector Umin = (Ui.size() ? U0 : vector::uniform(GREAT));
            vector Umax = (Ui.size() ? U0 : vector::uniform(-GREAT));
            reduce(Umin, minOp<vector>());
            reduce(Umax, maxOp<vector>());
            const bool nonUniform = dev > 0 || mag(Umax - Umin) > 0;
            startup.startProbe(!potentialInit && nonUniform);
        }
        scalar betaGlobal = startup.beta(iter);
        const bool startupDone = startup.done(iter);
        if (iter > 1 && betaGlobal != startup.beta(iter - 1))
        {
            // The discretisation changes along the ramp: the Anderson
            // history refers to another operator
            flushHistory("startupRamp");
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

        // Local limiter memory: release step (no effect without memory)
        ptc.beginIteration();

        // Level 3: dU_P/(fLoc Uref) per cell of the accepted assembly
        scalarField locRatio;
        if (diag.active(3))
        {
            locRatio.resize(mesh.nCells(), Zero);
        }

        // V/dt of the accepted assembly (SFD filter step)
        scalarField rDTVacc;

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
                assembler.momentumDiag(),
                ls.Uref(),
                (diag.active(3) ? &locRatio : nullptr)
            );
            dfields.record(rDTV, cflF, beta);                   // C6
            assembler.assembleContinuity(rDTV);
            tAsm += ta.elapsedTime();
            if (sfd.active())
            {
                rDTVacc = rDTV;
            }

            if (diag.active(1))
            {
                // Local pseudo-time step dt_P = V_P/(V_P/dt_P) after the
                // local limit: min/median/max on this rank
                const doubleScalar t0 = diagnostics::clock();
                const scalarField& V = mesh.V();
                std::vector<doubleScalar> dts(std::size_t(rDTV.size()));
                forAll(rDTV, celli)
                {
                    // GUARD: rDeltaTV >= VSMALL by construction (5.4)
                    dts[std::size_t(celli)] =
                        V[celli]/max(rDTV[celli], VSMALL);
                }
                if (!dts.empty())
                {
                    const auto mm = std::minmax_element(dts.begin(), dts.end());
                    dtMin = *mm.first;
                    dtMax = *mm.second;
                    const auto mid = std::ptrdiff_t(dts.size()/2);
                    std::nth_element(dts.begin(), dts.begin() + mid, dts.end());
                    dtMed = dts[std::size_t(mid)];
                }
                tDiagMain += diagnostics::clock() - t0;
            }

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
                    dfields.write(phi, sentinel::lastValidName(validIter));
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

            if (diag.active(1))
            {
                const doubleScalar t0 = diagnostics::clock();
                trialsJ.beginObject();
                trialsJ.add("CFL", doubleScalar(ptc.CFL()));
                trialsJ.add("eta", doubleScalar(eta));
                trialsJ.add("linIts", perf.nIterations);
                trialsJ.add("linConverged", perf.converged);
                trialsJ.add("linFinal", doubleScalar(perf.finalResidual));
                trialsJ.add("failed", solveFailed);
                trialsJ.add("omega", doubleScalar(omega));
                if (diag.active(3) && !solveFailed)
                {
                    // Physicality violations (rank-local) at the full step
                    // and at the line-search omega
                    label nU = 0, np = 0;
                    ls.countViolations(dx, 1, nU, np);
                    trialsJ.add("violU", nU);
                    trialsJ.add("violP", np);
                    ls.countViolations(dx, omega, nU, np);
                    trialsJ.add("violUAtOmega", nU);
                    trialsJ.add("violPAtOmega", np);
                }
                trialsJ.endObject();
                tDiagMain += diagnostics::clock() - t0;
            }

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
                flushHistory("omegaMin");
            }
            break;
        }

        // Anderson history is invalid after a CFL change (B5)
        if (cuts > 0 || skipStep)
        {
            flushHistory(skipStep ? "skipStep" : "cflCut");
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
        // C2: pseudo-inverse fallbacks of the tensorial Rhie-Chow D
        if (iter % coupledDefaults::rhieChowWarnInterval == 0)
        {
            const label nPinv = assembler.rc().takePseudoInverseWindow();
            if (nPinv)
            {
                WarningInFunction
                    << nPinv << " cells with a singular momentum block used"
                    << " the pseudo-inverse for the Rhie-Chow D in the last "
                    << coupledDefaults::rhieChowWarnInterval
                    << " iterations (run total "
                    << assembler.rc().nPseudoInverseTotal() << ", C2)"
                    << endl;
            }
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
                ++nSenChecks;
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
        if (diag.active(1))
        {
            const doubleScalar t0 = diagnostics::clock();
            assembler.rc().updateFlux(phi, U, p, assembler.noc(), false);
            MRF.makeRelative(phi);
            tFlux = diagnostics::clock() - t0;

            // Cell mass error sum_f phi_f of the new flux (rank-local)
            const doubleScalar t1 = diagnostics::clock();
            scalarField div(mesh.nCells(), Zero);
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            const scalarField& phiI = phi.primitiveField();
            forAll(own, facei)
            {
                div[own[facei]] += phiI[facei];
                div[nei[facei]] -= phiI[facei];
            }
            forAll(phi.boundaryField(), patchi)
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const scalarField& pp = phi.boundaryField()[patchi];
                forAll(fc, pf)
                {
                    div[fc[pf]] += pp[pf];
                }
            }
            for (const scalar d : div)
            {
                massErrMax = max(massErrMax, doubleScalar(mag(d)));
                massErrSum += mag(d);
            }
            tDiagMain += diagnostics::clock() - t1;
        }
        else
        {
            assembler.rc().updateFlux(phi, U, p, assembler.noc(), false);
            MRF.makeRelative(phi);
        }

        // --- Sentinel (9.3)
        bool rolledBack = false;
        auto rollback = [&](const sentinel::checkResult& chk)
        {
            sen.restore(U, p, phi, kPtr, omegaPtr, nutPtr);
            flushHistory("rollback");
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
                dfields.write(phi, sentinel::lastValidName(validIter));
                FatalErrorInFunction
                    << sen.consecutive() << " consecutive rollbacks (limit "
                    << sen.maxRollbacks() << ") at iteration " << iter
                    << ". Last valid fields (iteration " << validIter
                    << ") written to " << sentinel::lastValidName(validIter)
                    << "." << exit(FatalError);
            }
        };

        {
            ++nSenChecks;
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
                flushHistory("dynamicSet");
            }

            // --- Turbulence, segregated (5.8)
            clockTime tt;
            laminarTransport.correct();
            turbulence->correct();
            turbSolved = true;
            if (diag.active(1))
            {
                // Cells the bounds below will modify (rank-local)
                const doubleScalar t0 = diagnostics::clock();
                if (kPtr)
                {
                    for (const scalar v : kPtr->primitiveField())
                    {
                        nBoundK += (v < kMin ? 1 : 0);
                    }
                }
                if (omegaPtr)
                {
                    for (const scalar v : omegaPtr->primitiveField())
                    {
                        nBoundOmega += (v < omegaMinBound ? 1 : 0);
                    }
                }
                tDiagMain += diagnostics::clock() - t0;
            }
            nNutCapped = applyBounds();
            tTurb = tt.elapsedTime();

            ++nSenChecks;
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

                // SFD filter step and deactivation (7.6)
                sfd.update(U, rDTVacc);
                sfd.checkOff(R, iter);

                // Line-search beta (7.2): CFL boost on a full step
                if (omega >= 1 && cuts == 0)
                {
                    ptc.boost(ls.beta());
                }
            }

            if (startup.update(iter, R))
            {
                // beta changes from the next iteration on: the Anderson
                // history refers to the upwind operator
                flushHistory("startupEnd");
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
                flushHistory("autoTune");
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
        {
            // Smoother pivot fallbacks of the last solve (D-049)
            const blockGAMGPrecon* gpf =
                dynamic_cast<const blockGAMGPrecon*>
                (
                    linSolver->preconditioner()
                );
            if (gpf)
            {
                const label nFb =
                    returnReduce(gpf->gamg().nPivotFallback(), sumOp<label>());
                nPivotFallbackTotal += nFb;
                if (nFb)
                {
                    Info<< " nPivFb=" << nFb;
                }
            }
        }
        if (conv.haveForces())
        {
            Info<< " Cd=" << conv.Cd() << " Cl=" << conv.Cl();
            // Window mean and RMS fluctuation (D-045 b), appended at the end
            // so that existing parsers of the CF| line keep working
            forceHistIter.append(iter);
            for (label ci = 0; ci < conv.nCoeffs(); ++ci)
            {
                const convergenceMonitor::coeffStats st = conv.stats(ci);
                const char* cn = convergenceMonitor::coeffName(ci);
                Info<< ' ' << cn << "Mean=" << st.mean
                    << ' ' << cn << "Rms=" << st.rms;
                forceHistMean[ci].append(st.mean);
                forceHistRms[ci].append(st.rms);
                forceHistDrift[ci].append(st.drift);
            }
        }
        Info<< endl;

        // --- Diagnostics record (TASK 5), no collectives at any level here
        if (diag.active(1))
        {
            const doubleScalar tRec0 = diagnostics::clock();
            trialsJ.endArray();

            diagJson j;
            j.beginObject();
            j.add("type", "iter");
            j.add("iter", iter);
            j.add
            (
                "phase",
                phase.classify
                (
                    betaGlobal, R, conv.residualTol(),
                    conv.forceCriterionRatio()
                )
            );
            j.add("wallTime", doubleScalar(runTimer.elapsedTime()));

            j.beginObject("residuals");
            j.add("R", doubleScalar(R));
            j.add("Rraw", doubleScalar(Rraw));
            j.add("rU", doubleScalar(assembler.rU()));
            j.add("rp", doubleScalar(assembler.rp()));
            j.add("R1", doubleScalar(R1));
            j.add("normFactor", doubleScalar(assembler.normFactor()));
            j.add("linInitial", doubleScalar(perf.initialResidual));
            j.add("linFinal", doubleScalar(perf.finalResidual));
            j.add("linIts", perf.nIterations);
            j.add("linRestarts", perf.nRestarts);
            j.add("linConverged", perf.converged);
            j.add("linBreakdown", perf.breakdown);
            j.add("rho", doubleScalar(perf.rho));
            j.add("massErrMax", massErrMax);
            j.add("massErrSum", massErrSum);
            j.endObject();

            j.beginObject("controls");
            j.add("CFL", doubleScalar(ptc.CFL()));
            j.add("CFLstart", doubleScalar(CFLstart));
            // GUARD: CFL >= CFLmin > 0
            j.add
            (
                "growth",
                doubleScalar(ptc.CFL()/max(CFLstart, VSMALL))
            );
            j.add("strategy", std::string(ptc.strategyName()));
            j.add("hold", ptc.holdRemaining());
            j.add("nLocLim", nLocLim);
            j.add("nLocThrottled", ptc.nLocalThrottled());
            j.add("nLocSticky", ptc.nLocalSticky());
            j.beginObject("dt");
            j.add("min", dtMin);
            j.add("median", dtMed);
            j.add("max", dtMax);
            j.endObject();
            j.add("eta", doubleScalar(eta));
            j.add("etaRaw", doubleScalar(ew.lastRaw()));
            j.add("etaClip", ew.lastClip());
            j.add("omega", doubleScalar(omega));
            j.add("cuts", cuts);
            j.add("skipStep", skipStep);
            j.addRaw("trials", trialsJ.str());
            j.beginObject("sentinel");
            j.add("checks", nSenChecks);
            j.add("rolledBack", rolledBack);
            j.add("nRollbacks", sen.nRollbacks());
            j.add("consecutive", sen.consecutive());
            j.endObject();
            j.beginObject("remediation");
            j.add("nStat", rem.nStatic());
            j.add("nDyn", rem.nDynamic());
            j.add("nDynSticky", rem.nSticky());
            j.add("nDynRamping", rem.nRamping());
            j.add("version", rem.dynamicVersion());
            j.endObject();
            j.beginObject("sfd");
            j.add("enabled", sfd.enabled());
            j.add("active", sfd.active());
            j.add("chiStar", doubleScalar(sfd.chiStar()));
            j.add("resets", sfd.nResets());
            j.add
            (
                "maxDev",
                doubleScalar(sfd.maxDeviation(U, ls.Uref()))
            );
            j.endObject();
            j.beginObject("anderson");
            j.add("enabled", aa.enabled());
            if (aa.enabled())
            {
                doubleScalar g2 = 0;
                for (const doubleScalar g : aa.lastGamma())
                {
                    g2 += g*g;
                }
                j.add("status", andersonStatus);
                j.add("m", aa.m());
                j.add("nHistory", aa.nHistory());
                j.add("gammaNorm", std::sqrt(g2));
                j.add("maxAbsGamma", aa.lastMaxAlpha());
                const label nFl = aa.nFlushed() - nFlushedStart;
                j.add("flushes", nFl);
                if (nFl > 0)
                {
                    j.add("flushReason", aa.lastFlushReason());
                }
                else
                {
                    j.addNull("flushReason");
                }
            }
            j.endObject();
            j.add("beta", doubleScalar(betaGlobal));
            j.add("startupDone", startupDone);
            j.endObject();

            j.beginObject("turbulence");
            for (const char* fld : {"k", "omega"})
            {
                Pair<SolverPerformance<scalar>> sp;
                if
                (
                    turbSolved
                 && mesh.data().solverPerformanceDict().readIfPresent(word(fld), sp)
                )
                {
                    j.beginObject(fld);
                    j.add("init", doubleScalar(sp.first().initialResidual()));
                    j.add("final", doubleScalar(sp.second().finalResidual()));
                    j.add("its", label(sp.second().nIterations()));
                    j.endObject();
                }
                else
                {
                    j.addNull(fld);
                }
            }
            j.add("nBoundK", nBoundK);
            j.add("nBoundOmega", nBoundOmega);
            j.add("nNutCapped", nNutCapped);
            j.add("nClamped", nClamped);
            j.endObject();

            if (conv.haveForces())
            {
                j.beginObject("forces");
                j.add("window", conv.rmsWindow());
                j.add("driftEnabled", conv.driftEnabled());
                j.add("driftConverged", conv.driftConverged());
                for (label ci = 0; ci < conv.nCoeffs(); ++ci)
                {
                    const convergenceMonitor::coeffStats st = conv.stats(ci);
                    j.beginObject(convergenceMonitor::coeffName(ci));
                    j.add("last", doubleScalar(conv.last(ci)));
                    j.add("mean", doubleScalar(st.mean));
                    j.add("rms", doubleScalar(st.rms));
                    j.add("drift", doubleScalar(st.drift));
                    j.add("n", st.n);
                    j.endObject();
                }
                j.endObject();
            }
            else
            {
                j.addNull("forces");
            }

            {
                const coupledAssembler::timings& at = assembler.times();
                const doubleScalar tKrylov =
                    tSolve - diag.tPrecSetup() - diag.tPrecApply();
                j.beginObject("timings");
                j.add("tAsm", doubleScalar(tAsm));
                j.add("tMomentumOps", at.momentumOps);
                j.add("tBoundary", at.boundary);
                j.add("tContinuity", at.continuity);
                j.add("tRhieChow", at.rhieChow);
                j.add("tFlux", tFlux);
                j.add("tSolve", doubleScalar(tSolve));
                j.add("tPrecSetup", diag.tPrecSetup());
                j.add("tPrecApply", diag.tPrecApply());
                j.add("nPrecSetup", diag.nPrecSetup());
                j.add("nPrecApply", diag.nPrecApply());
                j.add("tKrylov", tKrylov);
                j.add("tTurb", doubleScalar(tTurb));
                // tDiag: diagnostics work of this iteration so far (hooks
                // in the solver/GAMG, the evaluations here, building this
                // record) plus the serialisation/write of the previous
                // record (tDiagWritePrev)
                j.add
                (
                    "tDiag",
                    diag.tDiag() + tDiagMain + diagWritePrev
                  + (diagnostics::clock() - tRec0)
                );
                j.add("tDiagWritePrev", diagWritePrev);
                j.add("tIter", doubleScalar(iterTimer.elapsedTime()));
                j.add("tWall", doubleScalar(runTime.elapsedClockTime()));
                j.endObject();
            }

            j.beginObject("gamg");
            if (gpDiag)
            {
                const blockGAMG& g = gpDiag->gamg();
                const bool changed = (g.hierarchyVersion() != diagHierVersion);
                j.add("nLevels", g.nLevels());
                j.add("Cop", doubleScalar(g.operatorComplexity()));
                j.add("cycle", std::string(blockGAMG::cycleName(g.cycleType())));
                j.add("nPostSweeps", g.nPostSweeps());
                j.add("setups", diag.nPrecSetup());
                j.add("reagglomerated", changed && diagHierVersion >= 0);
                j.add("hierarchyVersion", g.hierarchyVersion());
                if (changed)
                {
                    j.beginObject("hierarchy");
                    j.add("nLevels", g.nLevels());
                    j.addList("cellsPerLevel", g.globalCellsPerLevel());
                    j.addList("ranksPerLevel", g.ranksPerLevel());
                    List<doubleScalar> ratios(g.ratios().size());
                    forAll(ratios, i)
                    {
                        ratios[i] = g.ratios()[i];
                    }
                    j.addList("ratios", ratios);
                    j.add("Cop", doubleScalar(g.operatorComplexity()));
                    j.add("mergeLevels", g.mergeLevels());
                    j.add("denseCoarsest", g.denseCoarsest());
                    j.endObject();
                    diagHierVersion = g.hierarchyVersion();
                }
                else
                {
                    j.add("hierarchy", "unchanged");
                }
            }
            j.endObject();

            j.beginObject("memory");
            {
                // getrusage updates the peak lazily: never below current
                const label rss = diagnostics::currentRSSkB();
                j.add("rssKB", rss);
                j.add("peakRssKB", max(rss, diagnostics::peakRSSkB()));
            }
            j.endObject();

            j.endObject();

            // --- Level 3 items (rank-local)
            if (diag.active(3))
            {
                diagJson& e = diag.ext();

                // f. inflow/outflow switches of mixed-type patch faces
                e.beginObject("bcFlips");
                forAll(bcPatches, i)
                {
                    const label patchi = bcPatches[i];
                    const scalarField& pp = phi.boundaryField()[patchi];
                    boolList& prev = bcOut[i];
                    label nFlip = 0, nOut = 0;
                    if (!bcInit)
                    {
                        prev.resize(pp.size());
                    }
                    forAll(pp, pf)
                    {
                        const bool out = (pp[pf] > 0);
                        nOut += (out ? 1 : 0);
                        if (bcInit && out != prev[pf])
                        {
                            ++nFlip;
                        }
                        prev[pf] = out;
                    }
                    e.beginObject(mesh.boundary()[patchi].name().c_str());
                    e.add("nFaces", label(pp.size()));
                    e.add("nOutflow", nOut);
                    if (bcInit)
                    {
                        e.add("nFlips", nFlip);
                    }
                    else
                    {
                        e.addNull("nFlips");
                    }
                    e.endObject();
                }
                bcInit = true;
                e.endObject();

                // h. Locally CFL-limited cells (D-055): the diagTopLimited
                // cells with the largest dU_P/(fLoc Uref) of the accepted
                // assembly, their centres and memory factors, plus the
                // centroid of all limited cells (rank-local)
                {
                    const label nTop = coupledDefaults::diagTopLimited;
                    DynamicList<label> lim;
                    vector centroid(Zero);
                    forAll(locRatio, celli)
                    {
                        if (locRatio[celli] > 1)
                        {
                            lim.append(celli);
                            centroid += mesh.C()[celli];
                        }
                    }
                    if (lim.size())
                    {
                        centroid /= scalar(lim.size());
                    }
                    std::partial_sort
                    (
                        lim.begin(),
                        lim.begin() + min(nTop, lim.size()),
                        lim.end(),
                        [&](const label a, const label b)
                        {
                            return locRatio[a] > locRatio[b];
                        }
                    );
                    const scalarField& fLocal = ptc.localFactor();
                    const labelList& cLocal = ptc.localCount();
                    e.beginObject("localLimit");
                    e.add("nLimited", label(lim.size()));
                    e.beginArray("centroid");
                    for (direction c = 0; c < vector::nComponents; ++c)
                    {
                        e.value(doubleScalar(centroid[c]));
                    }
                    e.endArray();
                    e.beginArray("top");
                    for (label i = 0; i < min(nTop, lim.size()); ++i)
                    {
                        const label celli = lim[i];
                        e.beginObject();
                        e.add("cell", celli);
                        e.beginArray("C");
                        for (direction c = 0; c < vector::nComponents; ++c)
                        {
                            e.value(doubleScalar(mesh.C()[celli][c]));
                        }
                        e.endArray();
                        e.add("ratio", doubleScalar(locRatio[celli]));
                        if (fLocal.size())
                        {
                            e.add("f", doubleScalar(fLocal[celli]));
                            e.add("count", cLocal[celli]);
                        }
                        e.endObject();
                    }
                    e.endArray();
                    e.endObject();
                }

                // i. Dynamic-set members (D-055): the first diagTopLimited
                // cells of the set with centre, age and entry count
                {
                    const label nTop = coupledDefaults::diagTopLimited;
                    const labelList& age = rem.age();
                    const labelList& entries = rem.entries();
                    label n = 0;
                    e.beginObject("dynamicSet");
                    e.beginArray("cells");
                    forAll(age, celli)
                    {
                        if (age[celli] < 0)
                        {
                            continue;
                        }
                        if (n < nTop)
                        {
                            e.beginObject();
                            e.add("cell", celli);
                            e.beginArray("C");
                            for (direction c = 0; c < vector::nComponents; ++c)
                            {
                                e.value(doubleScalar(mesh.C()[celli][c]));
                            }
                            e.endArray();
                            e.add("age", age[celli]);
                            e.add("entries", entries[celli]);
                            e.endObject();
                        }
                        ++n;
                    }
                    e.endArray();
                    e.add("n", n);
                    e.endObject();
                }

                // g. Anderson internals
                if (aa.enabled())
                {
                    e.beginObject("andersonInternals");
                    e.add("conditionEstimate", aa.conditionEstimate());
                    List<doubleScalar> ag(aa.lastGamma().size());
                    forAll(ag, i)
                    {
                        ag[i] = std::fabs(aa.lastGamma()[i]);
                    }
                    e.addList("absGamma", ag);
                    e.add("maxAlpha", aa.maxAlpha());
                    e.add
                    (
                        "maxAlphaClip",
                        andersonStatus == label(anderson::status::skipped)
                     && aa.lastMaxAlpha() > aa.maxAlpha()
                    );
                    e.add("nSkippedTotal", aa.nSkipped());
                    e.endObject();
                }
            }

            const doubleScalar tW0 = diagnostics::clock();
            diag.writeRecord(j.str());
            diagWritePrev = diagnostics::clock() - tW0;
        }

        // --- Write
        if (!rolledBack && conv.converged(R))
        {
            converged = true;
            Info<< "coupledFoam: converged at iteration " << iter
                << " (R " << R << ")"
                << (conv.driftConverged() ? ", force-coefficient drift rule" : "")
                << endl;
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
        j.add("Uref", ls.Uref());
        j.add("pref", ls.pref());
        j.add("UrefMode", ls.UrefMode());
        j.add("UrefSource", ls.UrefSource());
        j.add("iterations", iter);
        j.add("converged", converged);
        j.add("finalR", lastR);
        j.add("finalCFL", ptc.CFL());
        j.add("lastLinearIterations", lastLinIters);
        j.add("cflCuts", nCflCutsTotal);
        j.add("linFailPolicy", linFailPolicy);
        j.add("linAcceptedUnconverged", nLinAccepted);
        j.add("pivotFallbacks", nPivotFallbackTotal);
        j.add("rollbacks", sen.nRollbacks());
        j.add("staticCells", rem.nStatic());
        j.add("dynamicCells", rem.nDynamic());
        j.add("dynamicStickyCells", rem.nSticky());
        j.add("localThrottledCells", ptc.nLocalThrottled());
        j.add("localStickyCells", ptc.nLocalSticky());
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
        j.add("rhieChowTensorial", assembler.rc().tensorial());
        j.add("nPseudoInverse", assembler.rc().nPseudoInverseTotal());
        j.add("ftz", ftzApplied);
        j.add("fpeTraps", runInfo::fpeActive());
        if (conv.haveForces())
        {
            j.add("Cd", conv.Cd());
            j.add("Cl", conv.Cl());

            // Window statistics (D-045 b): final values and histories
            diagJson fs;
            fs.beginObject();
            fs.add("window", conv.rmsWindow());
            fs.add("driftEnabled", conv.driftEnabled());
            fs.add("driftConverged", conv.driftConverged());
            for (label ci = 0; ci < conv.nCoeffs(); ++ci)
            {
                const convergenceMonitor::coeffStats st = conv.stats(ci);
                fs.beginObject(convergenceMonitor::coeffName(ci));
                fs.add("last", doubleScalar(conv.last(ci)));
                fs.add("mean", doubleScalar(st.mean));
                fs.add("rms", doubleScalar(st.rms));
                fs.add("drift", doubleScalar(st.drift));
                fs.add("n", st.n);
                fs.endObject();
            }
            fs.endObject();
            j.addRaw("forceStats", fs.str());

            diagJson fh;
            fh.beginObject();
            fh.addList("iter", forceHistIter);
            for (label ci = 0; ci < conv.nCoeffs(); ++ci)
            {
                const std::string cn(convergenceMonitor::coeffName(ci));
                fh.addList((cn + "Mean").c_str(), forceHistMean[ci]);
                fh.addList((cn + "Rms").c_str(), forceHistRms[ci]);
                fh.addList((cn + "Drift").c_str(), forceHistDrift[ci]);
            }
            fh.endObject();
            j.addRaw("forceHistory", fh.str());
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
