/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockBiCGStab.H"
#include "blockPreconditioner.H"
#include "doubleReduce.H"
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

    const label comm = matrix_.comm();
    const label n = matrix_.nRows();

    updatePreconditioner();

    const reduceScalar nf =
        (normFactorIn > 0 ? normFactorIn : normFactor(x, b));

    blockScalarList r(n);
    blockScalarList rA0(n);
    blockScalarList p(n, Zero);
    blockScalarList v(n, Zero);
    blockScalarList y(n);
    blockScalarList s(n);
    blockScalarList z(n);
    blockScalarList t(n);

    matrix_.residual(r, x, b);

    // GUARD: nf >= SMALL by construction (normFactor adds SMALL)
    perf.initialResidual = doubleReduce::norm2(r, comm)/nf;
    perf.finalResidual = perf.initialResidual;

    if (converged(perf.initialResidual, perf.initialResidual, 0))
    {
        perf.converged = true;
        return perf;
    }

    // Start (and restart) state
    auto restart = [&]()
    {
        rA0 = r;
        p = Zero;
        v = Zero;
    };
    restart();

    reduceScalar r0Sqr = doubleReduce::sumSqr(rA0, comm);
    reduceScalar rho = 1;
    reduceScalar alpha = 1;
    reduceScalar omega = 1;
    bool fresh = true;

    blockScalar* __restrict__ xp = x.data();
    blockScalar* __restrict__ rp = r.data();
    blockScalar* __restrict__ pp = p.data();
    blockScalar* __restrict__ vp = v.data();
    blockScalar* __restrict__ yp = y.data();
    blockScalar* __restrict__ sp = s.data();
    blockScalar* __restrict__ zp = z.data();
    blockScalar* __restrict__ tp = t.data();

    while (perf.nIterations < maxIter_)
    {
        const reduceScalar rhoOld = rho;
        rho = doubleReduce::dot(rA0, r, comm);

        // GUARD: breakdown of rho (spec 6.2)
        if (std::abs(rho) < doubleScalarVSMALL*r0Sqr)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            matrix_.residual(r, x, b);
            restart();
            r0Sqr = doubleReduce::sumSqr(rA0, comm);
            rho = 1;
            alpha = 1;
            omega = 1;
            fresh = true;
            continue;
        }

        if (fresh)
        {
            for (label i = 0; i < n; ++i)
            {
                pp[i] = rp[i];
            }
            fresh = false;
        }
        else
        {
            // GUARD: rhoOld and omega are non-zero (checked when set)
            const blockScalar beta = narrow((rho/rhoOld)*(alpha/omega));
            const blockScalar om = narrow(omega);
            for (label i = 0; i < n; ++i)
            {
                pp[i] = rp[i] + beta*(pp[i] - om*vp[i]);
            }
        }

        precondition(y, p);
        matrix_.Amul(v, y);

        const reduceScalar rA0v = doubleReduce::dot(rA0, v, comm);

        // GUARD: denominator of alpha
        if (std::abs(rA0v) < doubleScalarVSMALL*r0Sqr)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            matrix_.residual(r, x, b);
            restart();
            r0Sqr = doubleReduce::sumSqr(rA0, comm);
            rho = 1;
            alpha = 1;
            omega = 1;
            fresh = true;
            continue;
        }

        alpha = rho/rA0v;  // GUARD: rA0v checked above
        const blockScalar al = narrow(alpha);

        for (label i = 0; i < n; ++i)
        {
            sp[i] = rp[i] - al*vp[i];
        }

        ++perf.nIterations;

        const reduceScalar sRes = doubleReduce::norm2(s, comm)/nf;
        if (converged(sRes, perf.initialResidual, perf.nIterations))
        {
            for (label i = 0; i < n; ++i)
            {
                xp[i] += al*yp[i];
            }
            perf.finalResidual = sRes;
            perf.converged = true;
            return perf;
        }

        precondition(z, s);
        matrix_.Amul(t, z);

        const FixedList<reduceScalar, 2> tt_ts =
            doubleReduce::dot2(t, t, t, s, comm);

        // GUARD: omega denominator
        omega = tt_ts[1]/std::max(tt_ts[0], doubleScalarVSMALL);

        const blockScalar om = narrow(omega);
        for (label i = 0; i < n; ++i)
        {
            xp[i] += al*yp[i] + om*zp[i];
            rp[i] = sp[i] - om*tp[i];
        }

        perf.finalResidual = doubleReduce::norm2(r, comm)/nf;

        if (converged(perf.finalResidual, perf.initialResidual, perf.nIterations))
        {
            perf.converged = true;
            return perf;
        }

        // GUARD: omega breakdown (spec 6.2)
        if (std::abs(omega) < doubleScalarVSMALL)
        {
            if (perf.nRestarts >= maxRestarts_)
            {
                perf.breakdown = true;
                break;
            }
            ++perf.nRestarts;
            matrix_.residual(r, x, b);
            restart();
            r0Sqr = doubleReduce::sumSqr(rA0, comm);
            rho = 1;
            alpha = 1;
            omega = 1;
            fresh = true;
        }
    }

    return perf;
}


// ************************************************************************* //
