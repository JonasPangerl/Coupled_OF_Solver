/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGMRES.H"
#include "blockPreconditioner.H"
#include "doubleReduce.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockGMRES, 0);
    addToRunTimeSelectionTable(blockSolver, blockGMRES, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGMRES::blockGMRES
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    blockSolver(matrix, dict),
    restart_
    (
        max
        (
            label(1),
            dict.getOrDefault<label>
            (
                "gmresRestart",
                coupledDefaults::gmresRestart
            )
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockSolverPerformance Foam::blockGMRES::solve
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
    const label m = restart_;

    updatePreconditioner();

    const reduceScalar nf =
        (normFactorIn > 0 ? normFactorIn : normFactor(x, b));

    blockScalarList r(n);
    blockScalarList w(n);
    List<blockScalarList> V(m + 1);
    List<blockScalarList> Z(m);
    for (label j = 0; j <= m; ++j)
    {
        V[j].resize(n);
    }
    for (label j = 0; j < m; ++j)
    {
        Z[j].resize(n);
    }

    // Hessenberg (m+1) x m, Givens rotations, rhs of the LSQ problem
    List<reduceScalarList> H(m + 1, reduceScalarList(m, Zero));
    reduceScalarList cs(m, Zero);
    reduceScalarList sn(m, Zero);
    reduceScalarList g(m + 1, Zero);
    reduceScalarList yv(m, Zero);

    matrix_.residual(r, x, b);
    reduceScalar beta = doubleReduce::norm2(r, comm);

    perf.initialResidual = beta/nf;  // GUARD: nf >= SMALL
    perf.finalResidual = perf.initialResidual;

    if (converged(perf.initialResidual, perf.initialResidual, 0))
    {
        perf.converged = true;
        return perf;
    }

    while (perf.nIterations < maxIter_)
    {
        // v_0 = r/beta
        {
            // GUARD: beta > 0 here (otherwise converged above / below)
            const blockScalar rb = narrow(1.0/std::max(beta, doubleScalarVSMALL));
            blockScalar* __restrict__ v0 = V[0].data();
            const blockScalar* __restrict__ rp = r.cdata();
            for (label i = 0; i < n; ++i)
            {
                v0[i] = rb*rp[i];
            }
        }

        g = Zero;
        g[0] = beta;

        label j = 0;
        bool done = false;

        for (; j < m && perf.nIterations < maxIter_; ++j)
        {
            precondition(Z[j], V[j]);
            matrix_.Amul(w, Z[j]);

            // Modified Gram-Schmidt
            for (label i = 0; i <= j; ++i)
            {
                const reduceScalar hij = doubleReduce::dot(w, V[i], comm);
                H[i][j] = hij;
                const blockScalar h = narrow(hij);
                blockScalar* __restrict__ wp = w.data();
                const blockScalar* __restrict__ vi = V[i].cdata();
                for (label k = 0; k < n; ++k)
                {
                    wp[k] -= h*vi[k];
                }
            }

            const reduceScalar hNext = doubleReduce::norm2(w, comm);
            H[j + 1][j] = hNext;

            {
                // GUARD: lucky breakdown gives hNext = 0; V[j+1] unused then
                const blockScalar rh =
                    narrow(1.0/std::max(hNext, doubleScalarVSMALL));
                blockScalar* __restrict__ vn = V[j + 1].data();
                const blockScalar* __restrict__ wp = w.cdata();
                for (label k = 0; k < n; ++k)
                {
                    vn[k] = rh*wp[k];
                }
            }

            // Apply previous rotations to column j
            for (label i = 0; i < j; ++i)
            {
                const reduceScalar t = cs[i]*H[i][j] + sn[i]*H[i + 1][j];
                H[i + 1][j] = -sn[i]*H[i][j] + cs[i]*H[i + 1][j];
                H[i][j] = t;
            }

            // New rotation
            {
                const reduceScalar a = H[j][j];
                const reduceScalar bb = H[j + 1][j];
                const reduceScalar rr = std::sqrt(a*a + bb*bb);
                // GUARD: rotation denominator
                const reduceScalar d = std::max(rr, doubleScalarVSMALL);
                cs[j] = a/d;
                sn[j] = bb/d;
                H[j][j] = cs[j]*a + sn[j]*bb;
                H[j + 1][j] = 0;
                g[j + 1] = -sn[j]*g[j];
                g[j] = cs[j]*g[j];
            }

            ++perf.nIterations;
            perf.finalResidual = std::abs(g[j + 1])/nf;

            if
            (
                converged
                (
                    perf.finalResidual,
                    perf.initialResidual,
                    perf.nIterations
                )
             || H[j][j] == 0
            )
            {
                ++j;
                done = true;
                break;
            }
        }

        // Back substitution H(0:j,0:j) y = g(0:j)
        for (label i = j - 1; i >= 0; --i)
        {
            reduceScalar s = g[i];
            for (label k = i + 1; k < j; ++k)
            {
                s -= H[i][k]*yv[k];
            }
            // GUARD: diagonal of the rotated Hessenberg
            const reduceScalar hii = H[i][i];
            yv[i] =
                s/(std::abs(hii) > doubleScalarVSMALL
                 ? hii : doubleScalarVSMALL);
        }

        // x += Z y
        for (label i = 0; i < j; ++i)
        {
            const blockScalar yi = narrow(yv[i]);
            blockScalar* __restrict__ xp = x.data();
            const blockScalar* __restrict__ zi = Z[i].cdata();
            for (label k = 0; k < n; ++k)
            {
                xp[k] += yi*zi[k];
            }
        }

        // True residual for the restart (and as the reported value)
        matrix_.residual(r, x, b);
        beta = doubleReduce::norm2(r, comm);
        perf.finalResidual = beta/nf;

        // Only the true residual decides: the Givens estimate |g_j+1| can be
        // optimistic in float arithmetic. A cycle that ended early ("done")
        // but is not truly converged simply restarts.
        if (converged(perf.finalResidual, perf.initialResidual, perf.nIterations))
        {
            perf.converged = true;
            break;
        }

        if (done && beta <= doubleScalarVSMALL)
        {
            // Exact solution reached, nothing left to reduce
            break;
        }

        ++perf.nRestarts;
    }

    return perf;
}


void Foam::blockGMRES::writeSettings(dictionary& dict) const
{
    blockSolver::writeSettings(dict);
    dict.add("gmresRestart", restart_);
}


// ************************************************************************* //
