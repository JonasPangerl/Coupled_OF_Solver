/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockFGMRES.H"
#include "blockPreconditioner.H"
#include "doubleReduce.H"
#include "blockKernels.H"
#include "alignedList.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockFGMRES, 0);
    addToRunTimeSelectionTable(blockSolver, blockFGMRES, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockFGMRES::blockFGMRES
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
                "restart",
                dict.getOrDefault<label>("gmresRestart", coupledDefaults::restart)
            )
        )
    ),
    stagnationRatio_
    (
        dict.getOrDefault<doubleScalar>
        (
            "stagnationRatio",
            coupledDefaults::stagnationRatio
        )
    )
{
    // Negated comparison also rejects non-finite input
    if (!(stagnationRatio_ >= 0 && stagnationRatio_ < 1))
    {
        FatalIOErrorInFunction(dict)
            << "stagnationRatio must be in [0, 1) (0 = off), got "
            << stagnationRatio_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockSolverPerformance Foam::blockFGMRES::solve
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
    const label m = restart_;

    updatePreconditioner();

    const reduceScalar nf =
        (normFactorIn > 0 ? normFactorIn : normFactor(x, b));

    // Work vectors 64-byte aligned (6.5.3), kept across solves (D-069 F12)
    r_.resize_nocopy(n);
    w_.resize_nocopy(n);
    if (V_.size() != m + 1)
    {
        V_.clear();
        V_.resize(m + 1);
        Z_.clear();
        Z_.resize(m);
    }
    alignedList<blockScalar>& r = r_;
    alignedList<blockScalar>& w = w_;
    List<alignedList<blockScalar>>& V = V_;
    List<alignedList<blockScalar>>& Z = Z_;
    List<const blockScalar*> zPtr(m);
    for (label j = 0; j <= m; ++j)
    {
        V[j].resize_nocopy(n);
    }
    for (label j = 0; j < m; ++j)
    {
        Z[j].resize_nocopy(n);
        zPtr[j] = Z[j].cdata();
    }

    // Hessenberg (m+1) x m, Givens rotations, rhs of the LSQ problem
    List<reduceScalarList> H(m + 1, reduceScalarList(m, Zero));
    reduceScalarList cs(m, Zero);
    reduceScalarList sn(m, Zero);
    reduceScalarList g(m + 1, Zero);
    reduceScalarList yv(m, Zero);

    // Iterate in double (iterative refinement, see class description)
    xd_.resize_nocopy(n);
    alignedList<reduceScalar>& xd = xd_;
    blockKernels::widen(n, x.cdata(), xd.data());

    matrix_.residualDouble(r, xd, b);
    reduceScalar beta = doubleReduce::norm2(r, comm);
    const reduceScalar beta0 = beta;

    perf.initialResidual = beta/nf;  // GUARD: nf >= cfVSmall
    perf.finalResidual = perf.initialResidual;

    if (converged(perf.initialResidual, perf.initialResidual, 0))
    {
        perf.converged = true;
        return perf;
    }

    // True residual at the end of the previous restart cycle (D-080)
    doubleScalar cyclePrev = perf.initialResidual;

    while (perf.nIterations < maxIter_)
    {
        // v_0 = r/beta
        {
            // GUARD: beta > 0 here (otherwise converged above / below)
            const blockScalar rb = narrow(1.0/std::max(beta, doubleScalarVSMALL));
            blockKernels::scale(n, rb, r.cdata(), V[0].data());
        }

        g = Zero;
        g[0] = beta;

        label j = 0;
        bool done = false;

        for (; j < m && perf.nIterations < maxIter_; ++j)
        {
            precondition(Z[j], V[j]);
            matrix_.Amul(w, Z[j]);

            // rho on the first preconditioner application (6.3.5):
            // ||v0 - w||^2 = 1 - 2 v0.w + w.w  (v0 has unit norm); w.w is
            // taken in the same pass as h_0j
            const bool measure = (j == 0 && perf.rho < 0);
            reduceScalar ww = 0;

            // Modified Gram-Schmidt, fused (6.5.2): the pass that subtracts
            // h_ij v_i from w also forms h_(i+1)j = <w, v_(i+1)>, the last
            // pass forms ||w||^2. Same values as the unfused MGS (the dot
            // uses the updated w), one reduction per basis vector.
            reduceScalar hij = 0;
            if (measure)
            {
                const blockKernels::sumPair s =
                    blockKernels::dot_sumSqr(n, w.cdata(), V[0].cdata());
                reduceScalar sums[2] = {s.first, s.second};
                doubleReduce::parSum(sums, 2, comm);
                hij = sums[0];
                ww = sums[1];
            }
            else
            {
                hij = doubleReduce::parSum
                (
                    blockKernels::dot(n, w.cdata(), V[0].cdata()),
                    comm
                );
            }

            reduceScalar wwNext = 0;
            for (label i = 0; i <= j; ++i)
            {
                H[i][j] = hij;
                const blockScalar h = narrow(hij);
                if (i < j)
                {
                    hij = doubleReduce::parSum
                    (
                        blockKernels::mgs_axpy_dot
                        (
                            n, h, V[i].cdata(), w.data(), V[i + 1].cdata()
                        ),
                        comm
                    );
                }
                else
                {
                    wwNext = doubleReduce::parSum
                    (
                        blockKernels::update_residual_norm
                        (
                            n, h, V[i].cdata(), w.data()
                        ),
                        comm
                    );
                }
            }

            if (measure)
            {
                perf.rho = std::sqrt(max(1 - 2*H[0][0] + ww, 0.0));
                if (diagActive(2))
                {
                    // Scale-free rhoOpt = min_a ||v0 - a w|| (D-039), from
                    // the values of this pass: no extra reduction
                    // GUARD: w.w > 0 unless the preconditioned vector is 0
                    diag_->rhoOpt
                    (
                        std::sqrt
                        (
                            max
                            (
                                1 - H[0][0]*H[0][0]
                               /std::max(ww, doubleScalarVSMALL),
                                0.0
                            )
                        )
                    );
                }
            }

            // GUARD: sum of squares >= 0 by construction; clamp for sqrt
            const reduceScalar hNext = std::sqrt(wwNext > 0 ? wwNext : 0);
            H[j + 1][j] = hNext;

            // Happy breakdown (B1): the Krylov space is invariant
            const bool happy = hNext < doubleScalarVSMALL*beta0;

            {
                // GUARD: lucky breakdown gives hNext = 0; V[j+1] unused then
                const blockScalar rh =
                    narrow(1.0/std::max(hNext, doubleScalarVSMALL));
                blockKernels::scale(n, rh, w.cdata(), V[j + 1].data());
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
            if (diagActive(2))
            {
                diag_->krylovResidual(perf.finalResidual);
            }

            if (debug)
            {
                Info<< typeName << ": iteration " << perf.nIterations
                    << " Arnoldi estimate " << perf.finalResidual << endl;
            }

            if
            (
                converged
                (
                    perf.finalResidual,
                    perf.initialResidual,
                    perf.nIterations
                )
             || happy
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

        // x += Z y, in double (the products of a double and a blockScalar
        // value are formed in double); one chunked pass over x (6.5.2),
        // per entry in the order i = 0 .. j-1 as the former j passes
        blockKernels::axpyMultiDouble
        (
            n, j, yv.cdata(), zPtr.cdata(), xd.data()
        );

        // True residual of the double iterate, evaluated in double: seeds the
        // restart and is the reported value
        matrix_.residualDouble(r, xd, b);
        beta = doubleReduce::norm2(r, comm);
        perf.finalResidual = beta/nf;
        if (diagActive(2))
        {
            diag_->restartResidual(perf.finalResidual);
        }

        if (debug)
        {
            Info<< typeName << ": cycle end after " << j << " steps, true "
                << "residual " << perf.finalResidual << endl;
        }

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

        // D-080: a full cycle that did not reduce the true residual below
        // stagnationRatio x the previous cycle's will not converge within
        // maxIter either; stop it now, unconverged, so that the caller cuts
        // the CFL after a few cycles instead of after maxIter iterations
        if
        (
            stagnationRatio_ > 0
         && perf.finalResidual > stagnationRatio_*cyclePrev
        )
        {
            if (debug)
            {
                Info<< typeName << ": stagnation after " << perf.nIterations
                    << " iterations (" << perf.finalResidual << " > "
                    << stagnationRatio_ << " x " << cyclePrev << ")" << endl;
            }
            break;
        }
        cyclePrev = perf.finalResidual;

        ++perf.nRestarts;
    }

    // Return the double iterate rounded to blockScalar
    blockKernels::narrowCopy(n, xd.cdata(), x.data());

    return perf;
}


void Foam::blockFGMRES::writeSettings(dictionary& dict) const
{
    blockSolver::writeSettings(dict);
    dict.add("restart", restart_);
    dict.add("stagnationRatio", stagnationRatio_);
}


// ************************************************************************* //
