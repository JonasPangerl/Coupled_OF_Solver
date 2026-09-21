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
    const label startupUpwindIters =
        coupledDict.getOrDefault<label>
        (
            "startupUpwindIters",
            coupledDefaults::startupUpwindIters
        );
    const scalar startupSwitchR =
        coupledDict.getOrDefault<scalar>
        (
            "startupSwitchR",
            coupledDefaults::startupSwitchR
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
        eff.add("startupUpwindIters", startupUpwindIters);
        eff.add("startupSwitchR", startupSwitchR);
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
        eff.add("linearSolver", lin);

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
    bool startupDone = false;
    scalar R1 = -1;
    scalar phiConsistency = -1;
    bool restarted = false;

    {
        dictionary st;
        if (state.read(st))
        {
            restarted = true;
            iter = st.get<label>("iter");
            startupDone = st.get<bool>("startupDone");
            R1 = st.get<scalar>("R1");
            ls.setReference(st.get<scalar>("Uref"), st.get<scalar>("pref"));
            ptc.readState(st);
            rem.readState(st);
            sen.readState(st);

            Info<< "coupledFoam: restart from " << runTime.timeName()
                << " at iteration " << iter << ", CFL " << ptc.CFL()
                << ", startupDone " << startupDone << endl;

            // phi consistency with the stored D (spec 10)
            tmp<volScalarField> tD = state.readD();
            if (tD.valid())
            {
                const scalarField& Dc = tD().primitiveField();
                scalarField abar(mesh.V()/max(Dc, VSMALL));  // GUARD
                assembler.rc().updateD(abar, p);

                surfaceScalarField phiRe("phiRecomputed", phi);
                assembler.rc().updateFlux(phiRe, U, p, assembler.noc());
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
        st.set("startupDone", startupDone);
        st.set("R1", R1);
        st.set("Uref", ls.Uref());
        st.set("pref", ls.pref());
        ptc.writeState(st);
        rem.writeState(st);
        sen.writeState(st);
        st.set("refinementHistory", labelList());
        return st;
    };

    auto writeOutputs = [&]()
    {
        rem.write();
        state.writeD(assembler.rc().D());
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
        const scalarField beta(rem.beta(startupDone ? 1.0 : 0.0));
        const scalarField cflF(rem.cflFactor());

        blockScalarList dx(blockDim*mesh.nCells(), Zero);
        blockSolverPerformance perf;
        scalar omega = 1;
        label cuts = 0;
        label nLocLim = 0;
        scalar Rraw = 0;

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

            omega = ls.omega(dx);

            if (omega < ls.omegaMin() && cuts < ls.maxCflCuts())
            {
                ptc.decrease(ls.kappa());
                ++cuts;
                ++nCflCutsTotal;
                continue;
            }
            if (omega < ls.omegaMin())
            {
                omega = ls.omegaMin();
                rem.markDynamic(ls.offendingCells(dx));
            }
            break;
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

        // --- Field update
        rem.clipIncrement(dx, U, omega, ls.Uref());
        sen.store(U, p, phi, kPtr, omegaPtr, nutPtr);

        {
            vectorField& Ui = U.primitiveFieldRef();
            scalarField& pi = p.primitiveFieldRef();
            forAll(Ui, celli)
            {
                const blockScalar* d = dx.cdata() + celli*blockDim;
                Ui[celli] += omega*vector(d[0], d[1], d[2]);
                pi[celli] += omega*scalar(d[blockP]);
            }
        }
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        assembler.rc().updateFlux(phi, U, p, assembler.noc());
        MRF.makeRelative(phi);

        // --- Sentinel (9.3)
        bool rolledBack = false;
        auto rollback = [&](const sentinel::checkResult& chk)
        {
            sen.restore(U, p, phi, kPtr, omegaPtr, nutPtr);
            ptc.decrease(sen.cflFactor());
            rem.markDynamic(chk.offending);
            rolledBack = true;

            Info<< "coupledFoam: sentinel rollback at iteration " << iter
                << " (nonFinite " << chk.nNonFinite << ", max|U| "
                << chk.maxMagU << ", p [" << chk.minP << ", " << chk.maxP
                << "]), CFL -> " << ptc.CFL() << endl;

            if (sen.exhausted())
            {
                sen.writeLastValid
                (
                    U, p, phi, kPtr, omegaPtr, nutPtr, chk.offending, iter
                );
                rem.write();
                FatalErrorInFunction
                    << sen.consecutive() << " consecutive rollbacks (limit "
                    << sen.maxRollbacks() << ") at iteration " << iter
                    << ". Last valid fields written to " << iter
                    << "_lastValid." << exit(FatalError);
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
            rem.updateDynamic
            (
                U, p, kPtr, omegaPtr, ls.Uref(), ls.pref(), iter
            );

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
            sen.accepted();
            ptc.update(R);

            if
            (
                !startupDone
             && (iter >= startupUpwindIters || R < startupSwitchR)
            )
            {
                startupDone = true;
                Info<< "coupledFoam: start-up phase done at iteration "
                    << iter << " (R " << R << ")" << endl;
            }
        }

        conv.record(runTime);
        lastR = R;
        lastLinIters = perf.nIterations;

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
            << " nStat=" << rem.nStatic()
            << " nDyn=" << rem.nDynamic()
            << " nLocLim=" << nLocLim
            << " nRollback=" << sen.nRollbacks()
            << " nClamped=" << nClamped;
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
            runTime.writeNow();
            writeOutputs();
            break;
        }

        if (iter >= maxIter)
        {
            Info<< "coupledFoam: maxIter " << maxIter << " reached" << endl;
            runTime.writeNow();
            writeOutputs();
            break;
        }

        if (runTime.writeTime())
        {
            runTime.write();
            writeOutputs();
        }
    }

    // * * * * * * * * * * * * * * Run summary * * * * * * * * * * * * * * * //

    const label rssKB = runInfo::peakRSSkB();
    const label rssMax = returnReduce(rssKB, maxOp<label>());
    const scalar rssSum = returnReduce(scalar(rssKB), sumOp<scalar>());

    Info<< nl << "coupledFoam: iterations " << iter
        << ", converged " << converged
        << ", final R " << lastR
        << ", CFL cuts " << nCflCutsTotal
        << ", rollbacks " << sen.nRollbacks() << nl
        << "coupledFoam: wall time " << runTimer.elapsedTime() << " s" << nl
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
        j.add("rollbacks", sen.nRollbacks());
        j.add("staticCells", rem.nStatic());
        j.add("dynamicCells", rem.nDynamic());
        j.add("wallSeconds", runTimer.elapsedTime());
        j.add("peakRSS_MB_maxRank", rssMax/1024.0);
        j.add("peakRSS_MB_sum", rssSum/1024.0);
        j.add("restarted", restarted);
        j.add("phiConsistency", phiConsistency);
        j.add("ftz", ftzApplied);
        j.add("fpeTraps", runInfo::fpeActive());
        if (conv.haveForces())
        {
            j.add("Cd", conv.Cd());
            j.add("Cl", conv.Cl());
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
