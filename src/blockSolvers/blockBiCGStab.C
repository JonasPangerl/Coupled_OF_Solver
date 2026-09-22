/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockBiCGStab.H"
#include "blockPreconditioner.H"
#include "doubleReduce.H"
#include "blockKernels.H"
#include "alignedList.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockBiCGStab, 0);
    addToRunTimeSelectionTable(blockSolver, blockBiCGStab, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockBiCGStab::blockBiCGStab
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    blockSolver(matrix, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockSolverPerformance Foam::blockBiCGStab::solve
(
    blockScalarUList& x,
    const blockScalarUList& b,
    const reduceScalar normFactorIn
) const
{
    blockSolverPerformance perf;
    perf.solverName = typeName;
    const diagSolveScope diagScope(*this, perf);

    const label comm = matrix_.comm();
    const label n = matrix_.nRows();

    updatePreconditioner();

    const reduceScalar nf =
        (normFactorIn > 0 ? normFactorIn : normFactor(x, b));

    // Work vectors 64-byte aligned (6.5.3). s = r - alpha v is formed in
    // place in r (update_residual_norm), so no separate s vector.
    alignedList<blockScalar> r(n);
    alignedList<blockScalar> rA0(n);
    alignedList<blockScalar> p(n, Zero);
    alignedList<blockScalar> v(n, Zero);
    alignedList<blockScalar> y(n);
    alignedList<blockScalar> z(n);
    alignedList<blockScalar> t(n);

    matrix_.residual(r, x, b);

    // GUARD: nf >= cfVSmall by construction (normFactor adds cfVSmall)
    perf.initialResidual = doubleReduce::norm2(r, comm)/nf;
    perf.finalResidual = perf.initialResidual;

    if (converged(perf.initialResidual, perf.initialResidual, 0))
    {
        perf.converged = true;
        return perf;
    }

    // Best iterate: BiCGStab convergence is not monotone; near the float
    // floor it can diverge. The iterate with the smallest residual is kept
    // and returned if the solve does not converge.
    alignedList<blockScalar> xBest(x);
    reduceScalar bestRes = perf.initialResidual;

    auto keepBest = [&](const reduceScalar res)
    {
        if (res < bestRes)
        {
            bestRes = res;
            xBest.copyFrom(x);
        }
    };

    auto finish = [&]() -> blockSolverPerformance&
    {
        if (!perf.converged && bestRes < perf.finalResidual)
        {
            blockKernels::copy(n, xBest.cdata(), x.data());
            perf.finalResidual = bestRes;
        }
        return perf;
    };

    // Start (and restart) state
    auto restart = [&]()
    {
        rA0.copyFrom(r);
        p = Zero;
        v = Zero;
    };
    restart();

    reduceScalar r0Sqr = doubleReduce::sumSqr(rA0, comm);
    reduceScalar rho = 1;
    reduceScalar alpha = 1;
    reduceScalar omega = 1;
    bool fresh = true;

    // <rA0, r> of the current r, known from the fused update of r at the
    // end of the previous iteration (axpy_dot); recomputed after a restart
    reduceScalar rhoNext = 0;
    bool rhoKnown = false;

    blockScalar* __restrict__ xp = x.data();
    blockScalar* __restrict__ rp = r.data();
    blockScalar* __restrict__ pp = p.data();
    blockScalar* __restrict__ vp = v.data();
    const blockScalar* __restrict__ yp = y.cdata();
    const blockScalar* __restrict__ zp = z.cdata();
    const blockScalar* __restrict__ tp = t.cdata();
    const blockScalar* __restrict__ r0p = rA0.cdata();

    auto restartState = [&]()
    {
        matrix_.residual(r, x, b);
        restart();
        r0Sqr = doubleReduce::sumSqr(rA0, comm);
        rho = 1;
        alpha = 1;
        omega = 1;
        fresh = true;
        rhoKnown = false;
    };

    while (perf.nIterations < maxIter_)
    {
        const reduceScalar rhoOld = rho;
        rho =
            rhoKnown
          ? rhoNext
          : doubleReduce::parSum(blockKernels::dot(n, r0p, rp), comm);
        rhoKnown = false;

        // GUARD: breakdown of rho (spec 6.2)
        if (std::abs(rho) < doubleScalarVSMALL*r0Sqr)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            restartState();
            continue;
        }

        if (fresh)
        {
            blockKernels::copy(n, rp, pp);
            fresh = false;
        }
        else
        {
            // GUARD: rhoOld and omega are non-zero (checked when set)
            const blockScalar beta = narrow((rho/rhoOld)*(alpha/omega));
            const blockScalar om = narrow(omega);
            #pragma omp simd
            for (label i = 0; i < n; ++i)
            {
                pp[i] = rp[i] + beta*(pp[i] - om*vp[i]);
            }
        }

        precondition(y, p);
        matrix_.Amul(v, y);

        // rho on the first preconditioner application of the solve (6.3.5):
        // p = r0 here, rho = ||r0 - A M^-1 r0|| / ||r0||
        if (perf.rho < 0)
        {
            reduceScalar d2 = 0;
            reduceScalar p2 = 0;
            #pragma omp simd reduction(+:d2, p2)
            for (label i = 0; i < n; ++i)
            {
                const reduceScalar pi = toDouble(pp[i]);
                const reduceScalar di = pi - toDouble(vp[i]);
                d2 += di*di;
                p2 += pi*pi;
            }
            reduceScalar sums[2] = {d2, p2};
            doubleReduce::parSum(sums, 2, comm);
            // GUARD: ||r0|| > 0 (not converged at iteration 0)
            perf.rho =
                std::sqrt(sums[0])
               /max(std::sqrt(sums[1]), doubleScalarVSMALL);
        }

        const reduceScalar rA0v =
            doubleReduce::parSum(blockKernels::dot(n, r0p, vp), comm);

        // GUARD: denominator of alpha
        if (std::abs(rA0v) < doubleScalarVSMALL*r0Sqr)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            restartState();
            continue;
        }

        alpha = rho/rA0v;  // GUARD: rA0v checked above
        const blockScalar al = narrow(alpha);

        // s = r - alpha v in place, with ||s||^2 (update_residual_norm)
        const reduceScalar ss = doubleReduce::parSum
        (
            blockKernels::update_residual_norm(n, al, vp, rp),
            comm
        );

        ++perf.nIterations;

        // GUARD: sum of squares >= 0 by construction; clamp for sqrt
        const reduceScalar sRes = std::sqrt(ss > 0 ? ss : 0)/nf;
        if (converged(sRes, perf.initialResidual, perf.nIterations))
        {
            if (diagActive(2))
            {
                // One entry per iteration: the half-step residual only
                // when the iteration ends there
                diag_->krylovResidual(sRes);
            }
            #pragma omp simd
            for (label i = 0; i < n; ++i)
            {
                xp[i] += al*yp[i];
            }
            perf.finalResidual = sRes;
            perf.converged = true;
            return perf;
        }

        // r holds s from here on
        precondition(z, r);
        matrix_.Amul(t, z);

        // (<t,s>, <t,t>) in one pass
        const blockKernels::sumPair ts_tt =
            blockKernels::dot_sumSqr(n, tp, rp);
        reduceScalar tsums[2] = {ts_tt.first, ts_tt.second};
        doubleReduce::parSum(tsums, 2, comm);

        // GUARD: omega denominator
        omega = tsums[0]/std::max(tsums[1], doubleScalarVSMALL);

        const blockScalar om = narrow(omega);
        #pragma omp simd
        for (label i = 0; i < n; ++i)
        {
            xp[i] += al*yp[i] + om*zp[i];
        }

        // r = s - omega t, fused with ||r||^2 and <r, rA0> (axpy_dot); the
        // latter is rho of the next iteration
        const blockKernels::sumPair rr_rr0 =
            blockKernels::axpy_dot(n, -om, tp, rp, r0p);
        reduceScalar rsums[2] = {rr_rr0.first, rr_rr0.second};
        doubleReduce::parSum(rsums, 2, comm);
        rhoNext = rsums[1];
        rhoKnown = true;

        // GUARD: sum of squares >= 0 by construction; clamp for sqrt
        perf.finalResidual = std::sqrt(rsums[0] > 0 ? rsums[0] : 0)/nf;
        if (diagActive(2))
        {
            diag_->krylovResidual(perf.finalResidual);
        }

        if (converged(perf.finalResidual, perf.initialResidual, perf.nIterations))
        {
            perf.converged = true;
            return perf;
        }

        keepBest(perf.finalResidual);

        // GUARD: omega breakdown (spec 6.2)
        if (std::abs(omega) < doubleScalarVSMALL)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            restartState();
        }
    }

    return finish();
}


// ************************************************************************* //
