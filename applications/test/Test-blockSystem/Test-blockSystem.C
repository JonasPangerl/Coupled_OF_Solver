/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-blockSystem

Description
    Offline preconditioner study on a linear system written by coupledFoam
    (coupled.dumpLinearSystem (iterations), serial): diag, upper, lower,
    rhs, normFactor, eta and CFL of one outer iteration.

    Every variant of the -variants dictionary is merged (recursively) onto
    solvers.coupled of the case and solved twice with a zero initial guess:
      solve 1  relTol = the dumped eta (as in the run), includes the
               preconditioner set-up (agglomeration, coarse operators,
               smoother factorisation)
      solve 2  relTol = -deepRelTol (default 1e-6), preconditioner reused:
               asymptotic convergence
    Reported per variant: iterations, rho (unscaled, as coupledFoam),
    rhoOpt = min_a ||r - a A M^-1 r||/||r|| (scale-free), wall times.
    autoTune is switched off so that the configured cycle is measured.

Usage
    Test-blockSystem -system <file> -variants <dictFile> [-deepRelTol x]
        [-only name]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "blockSolver.H"
#include "blockPreconditioner.H"
#include "blockGAMGPrecon.H"
#include "doubleReduce.H"
#include "clockTime.H"
#include "IFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("system", "file", "Dumped linear system");
    argList::addOption("variants", "file", "Dictionary of variants");
    argList::addOption("deepRelTol", "value", "relTol of solve 2 (1e-6)");
    argList::addOption("only", "name", "Run only this variant");
    argList::addBoolOption("dominance", "Row dominance statistics");
    argList::addBoolOption("equilibrate", "Symmetric diagonal scaling");
    argList::addBoolOption("diagnose", "Worst cells of one application");
    argList::addOption("scalars", "(nf eta cfl)", "Override dumped scalars");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    blockLduMatrix4 A(mesh);
    blockScalarList rhs;
    doubleScalar normFactor = 0, eta = 0, cfl = 0;
    {
        IFstream is
        (
            args.get<fileName>("system"),
            IOstreamOption(IOstreamOption::BINARY)
        );
        // Header-less binary file of float lists: raw 4-byte values
        is.setScalarByteSize(sizeof(blockScalar));
        is  >> A.diag() >> A.upper() >> A.lower() >> rhs;
        if (args.found("scalars"))
        {
            // normFactor eta CFL given (dumps without separators)
            const List<doubleScalar> s(args.getList<doubleScalar>("scalars"));
            normFactor = s[0];
            eta = s[1];
            cfl = s[2];
        }
        else
        {
            is  >> normFactor >> eta >> cfl;
        }
    }
    if
    (
        A.diag().size() != blockSize*mesh.nCells()
     || A.upper().size() != blockSize*mesh.nInternalFaces()
     || rhs.size() != blockDim*mesh.nCells()
    )
    {
        FatalErrorInFunction
            << "Dumped system does not match the mesh" << exit(FatalError);
    }
    // Optional symmetric equilibration S A S y = S b, x = S y with
    // s = 1/sqrt|D_kk| per cell and component
    blockScalarList sc(rhs.size(), blockScalar(1));
    blockLduMatrix4 A0(mesh);
    A0.diag() = A.diag();
    A0.upper() = A.upper();
    A0.lower() = A.lower();
    A0.markUpdated();
    const blockScalarList rhs0(rhs);
    if (args.found("equilibrate"))
    {
        const label nc = mesh.nCells();
        for (label c = 0; c < nc; ++c)
        {
            for (label k = 0; k < blockDim; ++k)
            {
                const doubleScalar d =
                    std::abs(toDouble(A.diag()[c*blockSize + k*blockDim + k]));
                sc[c*blockDim + k] =
                    narrow(1.0/std::sqrt(std::max(d, doubleScalarVSMALL)));
            }
        }
        auto scaleBlock = [&](blockScalar* B, const label ri, const label ci)
        {
            for (label r = 0; r < blockDim; ++r)
            {
                for (label c = 0; c < blockDim; ++c)
                {
                    B[r*blockDim + c] *=
                        sc[ri*blockDim + r]*sc[ci*blockDim + c];
                }
            }
        };
        for (label c = 0; c < nc; ++c)
        {
            scaleBlock(A.diag().data() + c*blockSize, c, c);
        }
        const labelUList& l = mesh.lduAddr().lowerAddr();
        const labelUList& u = mesh.lduAddr().upperAddr();
        forAll(l, f)
        {
            scaleBlock(A.upper().data() + f*blockSize, l[f], u[f]);
            scaleBlock(A.lower().data() + f*blockSize, u[f], l[f]);
        }
        forAll(rhs, i)
        {
            rhs[i] *= sc[i];
        }
        normFactor = -1;
        Info<< "equilibrated" << nl;
    }
    A.markUpdated();

    // Original-norm reduction ||b - A0 S y||/||b|| of a scaled solution y
    auto trueRed = [&](const blockScalarList& y)
    {
        blockScalarList xx(y.size()), rr(y.size());
        forAll(y, i)
        {
            xx[i] = y[i]*sc[i];
        }
        A0.residual(rr, xx, rhs0);
        return doubleReduce::norm2(rr, A0.comm())
            /std::max(doubleReduce::norm2(rhs0, A0.comm()), doubleScalarVSMALL);
    };

    Info<< "system: " << args.get<fileName>("system") << " cells "
        << mesh.nCells() << " normFactor " << normFactor << " eta " << eta
        << " CFL " << cfl << nl << endl;

    if (args.found("dominance"))
    {
        // Row dominance: per row k, sum over neighbours of |A_kc| split into
        // u-columns (c<3) and the p-column, divided by |D_kk|, and the
        // in-block coupling sum_{c!=k}|D_kc|/|D_kk|
        const label nc = mesh.nCells();
        List<FixedList<doubleScalar, 4>> offU(nc), offP(nc);
        for (auto& v : offU) v = Zero;
        for (auto& v : offP) v = Zero;
        const labelUList& l = mesh.lduAddr().lowerAddr();
        const labelUList& u = mesh.lduAddr().upperAddr();
        auto acc = [&](const blockScalar* B, const label row)
        {
            for (label r = 0; r < blockDim; ++r)
            {
                for (label c = 0; c < blockDim; ++c)
                {
                    const doubleScalar v = std::abs(toDouble(B[r*blockDim + c]));
                    if (c < 3) offU[row][r] += v; else offP[row][r] += v;
                }
            }
        };
        forAll(l, f)
        {
            acc(A.upper().cdata() + f*blockSize, l[f]);
            acc(A.lower().cdata() + f*blockSize, u[f]);
        }
        for (label k = 0; k < blockDim; ++k)
        {
            doubleScalar mxU = 0, mxP = 0, mnD = GREAT, mxIn = 0;
            label nBad = 0;
            label cBad = -1;
            for (label c = 0; c < nc; ++c)
            {
                const doubleScalar d =
                    std::abs(toDouble(A.diag()[c*blockSize + k*blockDim + k]));
                doubleScalar in = 0;
                for (label j = 0; j < blockDim; ++j)
                {
                    if (j != k)
                    {
                        in += std::abs(toDouble(A.diag()[c*blockSize + k*blockDim + j]));
                    }
                }
                const doubleScalar rU = offU[c][k]/max(d, VSMALL);
                const doubleScalar rP = offP[c][k]/max(d, VSMALL);
                if (rU + rP > 1.0001) { ++nBad; }
                if (rU > mxU) { mxU = rU; cBad = c; }
                mxP = max(mxP, rP);
                mxIn = max(mxIn, in/max(d, VSMALL));
                mnD = min(mnD, d);
            }
            Info<< "DOM row " << k << " max offU/d " << mxU << " (cell "
                << cBad << " " << mesh.C()[cBad] << ") max offP/d " << mxP
                << " max inBlock/d " << mxIn << " min d " << mnD
                << " rows not dominant " << nBad << endl;
        }
    }

    const doubleScalar deepRelTol =
        args.getOrDefault<doubleScalar>("deepRelTol", 1e-6);

    const dictionary base(mesh.solversDict().subDict("coupled"));
    const IOdictionary variants
    (
        IOobject
        (
            args.get<fileName>("variants"),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        IFstream(args.get<fileName>("variants"))()
    );

    const word only(args.getOrDefault<word>("only", word::null));

    Info<< "RESULT variant nIter1 conv1 rho rhoOpt t1 nIter2 conv2 "
        << "res2 t2 tPerIt2" << endl;

    for (const entry& e : variants)
    {
        if (!e.isDict() || (!only.empty() && e.keyword() != only))
        {
            continue;
        }
        dictionary sd(base);
        sd.merge(e.dict());
        {
            dictionary g(sd.subOrEmptyDict("blockGAMG"));
            g.set("autoTune", false);
            sd.set("blockGAMG", g);
        }
        sd.set("maxIter", sd.getOrDefault<label>("maxIter", 200));

        autoPtr<blockSolver> solver = blockSolver::New(A, sd);

        // Solve 1: the run's relTol, includes the set-up
        solver->setRelTol(eta);
        blockScalarList x(rhs.size(), Zero);
        clockTime t1;
        const blockSolverPerformance p1 = solver->solve(x, rhs, normFactor);
        const doubleScalar s1 = t1.elapsedTime();
        const doubleScalar red1 = trueRed(x);

        // Scale-free first-application efficiency on r0 = rhs
        doubleScalar rhoOpt = -1;
        if (solver->preconditioner())
        {
            blockScalarList z(rhs.size(), Zero), w(rhs.size(), Zero);
            solver->preconditioner()->precondition(z, rhs);
            A.Amul(w, z);
            const doubleScalar rr = doubleReduce::sumSqr(rhs, A.comm());
            const doubleScalar ww = doubleReduce::sumSqr(w, A.comm());
            const doubleScalar rw = doubleReduce::dot(rhs, w, A.comm());
            rhoOpt = std::sqrt
            (
                std::max(1 - rw*rw/std::max(rr*ww, doubleScalarVSMALL), 0.0)
            );
        }

        if (solver->preconditioner())
        {
            // Timing: preconditioner set-up and one application
            blockPreconditioner& pc =
                const_cast<blockPreconditioner&>(*solver->preconditioner());
            const label nRep = 5;
            clockTime tu;
            for (label i = 0; i < nRep; ++i)
            {
                pc.update();
            }
            const doubleScalar tUpd = tu.elapsedTime()/nRep;
            blockScalarList z(rhs.size(), Zero);
            clockTime tp;
            for (label i = 0; i < nRep; ++i)
            {
                pc.precondition(z, rhs);
            }
            const doubleScalar tApp = tp.elapsedTime()/nRep;
            Info<< "TIMING " << e.keyword() << " update " << tUpd
                << " apply " << tApp << endl;
        }

        if (args.found("diagnose") && solver->preconditioner())
        {
            // Where does one preconditioner application go wrong?
            blockScalarList z(rhs.size(), Zero), w(rhs.size(), Zero);
            solver->preconditioner()->precondition(z, rhs);
            A.Amul(w, z);
            const label nc = mesh.nCells();
            List<doubleScalar> bad(nc, Zero);
            doubleScalar tot = 0, totR = 0;
            forAll(bad, c)
            {
                for (label k = 0; k < blockDim; ++k)
                {
                    const doubleScalar d =
                        toDouble(w[c*blockDim + k])
                      - toDouble(rhs[c*blockDim + k]);
                    bad[c] += d*d;
                    totR += sqr(toDouble(rhs[c*blockDim + k]));
                }
                tot += bad[c];
            }
            const labelList order(sortedOrder(bad));
            Info<< "DIAG ||Az-r||^2 " << tot << " ||r||^2 " << totR << nl;
            doubleScalar acc = 0;
            for (label i = 0; i < 10; ++i)
            {
                const label c = order[nc - 1 - i];
                acc += bad[c];
                Info<< "DIAG cell " << c << " C " << mesh.C()[c]
                    << " frac " << bad[c]/tot << " z (";
                for (label k = 0; k < blockDim; ++k)
                {
                    Info<< ' ' << z[c*blockDim + k];
                }
                Info<< ") r (";
                for (label k = 0; k < blockDim; ++k)
                {
                    Info<< ' ' << rhs[c*blockDim + k];
                }
                Info<< ") D (";
                for (label k = 0; k < blockSize; ++k)
                {
                    Info<< ' ' << A.diag()[c*blockSize + k];
                }
                Info<< ")" << nl;
            }
            Info<< "DIAG top10 fraction " << acc/tot << endl;
        }

        // Solve 2: deep, preconditioner reused
        solver->setRelTol(deepRelTol);
        x = Zero;
        clockTime t2;
        const blockSolverPerformance p2 = solver->solve(x, rhs, normFactor);
        const doubleScalar s2 = t2.elapsedTime();
        const doubleScalar red2 = trueRed(x);

        Info<< "RESULT " << e.keyword()
            << ' ' << p1.nIterations << ' ' << p1.converged
            << ' ' << p1.rho << ' ' << rhoOpt << ' ' << s1
            << ' ' << p2.nIterations << ' ' << p2.converged
            << ' ' << p2.finalResidual/max(p2.initialResidual, VSMALL)
            << ' ' << s2
            << ' ' << s2/max(p2.nIterations, label(1))
            << " trueRed1 " << red1 << " trueRed2 " << red2 << endl;
    }

    Info<< "End" << endl;
    return 0;
}


// ************************************************************************* //
