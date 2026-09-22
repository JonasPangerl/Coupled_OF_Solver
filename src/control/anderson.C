/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "anderson.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include <cmath>
#include <limits>

// * * * * * * * * * * * * * * * Local Constants * * * * * * * * * * * * * * //

namespace
{
    //- Unknowns per cell (u, v, w, p) and the index of p
    constexpr Foam::label nCmpt = Foam::vector::nComponents + 1;
    constexpr Foam::label pCmpt = Foam::vector::nComponents;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anderson::anderson(const fvMesh& mesh, const dictionary& coupledDict)
:
    mesh_(mesh),
    enabled_(coupledDefaults::andersonEnabled),
    m_(coupledDefaults::andersonM),
    beta_(coupledDefaults::andersonBeta),
    maxAlpha_(coupledDefaults::andersonMaxAlpha),
    n_(nCmpt*mesh.nCells()),
    nHist_(0),
    hasPrev_(false),
    sU_(0),
    sP_(0),
    nApplied_(0),
    nSkipped_(0),
    nRejected_(0),
    nFlushed_(0),
    lastStatus_(status::disabled),
    lastMaxAlpha_(0),
    lastGamma_(),
    lastFlushReason_("none")
{
    const dictionary& d = coupledDict.subOrEmptyDict("anderson");
    enabled_ = d.getOrDefault<bool>("enabled", coupledDefaults::andersonEnabled);
    m_ = d.getOrDefault<label>("m", coupledDefaults::andersonM);
    beta_ = d.getOrDefault<doubleScalar>("beta", coupledDefaults::andersonBeta);
    maxAlpha_ =
        d.getOrDefault<doubleScalar>
        (
            "maxAlpha",
            coupledDefaults::andersonMaxAlpha
        );

    if (m_ < 1 || !(beta_ > 0) || !(maxAlpha_ > 0))
    {
        FatalIOErrorInFunction(d)
            << "anderson: need m >= 1, beta > 0, maxAlpha > 0 (got m " << m_
            << ", beta " << beta_ << ", maxAlpha " << maxAlpha_ << ")"
            << exit(FatalIOError);
    }

    // B7 memory table: "Anderson m=4 (optional) +11.5 GB -> not enabled
    // above 35 M cells"
    if (enabled_)
    {
        const label nTotal = returnReduce(mesh.nCells(), sumOp<label>());
        if (nTotal > coupledDefaults::andersonMaxCells)
        {
            WarningInFunction
                << "anderson.enabled ignored: " << nTotal << " cells > "
                << coupledDefaults::andersonMaxCells
                << " (amendment B7 memory rule); Anderson acceleration is"
                << " disabled" << endl;
            enabled_ = false;
        }
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::anderson::sumAll(List<doubleScalar>& v)
{
    if (v.size())
    {
        reduce
        (
            v.data(),
            static_cast<int>(v.size()),
            sumOp<doubleScalar>(),
            UPstream::msgType(),
            UPstream::worldComm
        );
    }
}


void Foam::anderson::allocate()
{
    if (Q_.size() == m_ && xPrev_.size() == n_)
    {
        return;
    }

    Q_ = List<List<doubleScalar>>(m_, List<doubleScalar>(n_, Zero));
    D_ = List<List<doubleScalar>>(m_, List<doubleScalar>(n_, Zero));
    R_ = List<doubleScalar>(m_*m_, Zero);
    xPrev_ = List<doubleScalar>(n_, Zero);
    fPrev_ = List<doubleScalar>(n_, Zero);
    nHist_ = 0;
    hasPrev_ = false;
}


void Foam::anderson::removeOldest()
{
    const label n = nHist_;

    // Drop column 0 of R: R becomes upper Hessenberg n x (n-1)
    for (label j = 0; j < n - 1; ++j)
    {
        for (label i = 0; i <= j + 1; ++i)
        {
            R(i, j) = R(i, j + 1);
        }
    }
    for (label i = 0; i < n; ++i)
    {
        R(i, n - 1) = 0;
    }

    // Givens rotations restore the triangle; Q <- Q G^T
    for (label i = 0; i < n - 1; ++i)
    {
        const doubleScalar a = R(i, i);
        const doubleScalar b = R(i + 1, i);
        const doubleScalar r = std::hypot(a, b);

        // GUARD: r >= |b| = former diagonal R(i+2, i+2) > 0 (accepted
        // columns have a positive diagonal); skip if nothing to rotate
        if (!(r > 0))
        {
            continue;
        }
        const doubleScalar c = a/r;
        const doubleScalar s = b/r;

        for (label j = i; j < n - 1; ++j)
        {
            const doubleScalar ri = R(i, j);
            const doubleScalar ri1 = R(i + 1, j);
            R(i, j) = c*ri + s*ri1;
            R(i + 1, j) = -s*ri + c*ri1;
        }
        R(i + 1, i) = 0;

        List<doubleScalar>& qi = Q_[i];
        List<doubleScalar>& qi1 = Q_[i + 1];
        forAll(qi, k)
        {
            const doubleScalar a0 = qi[k];
            const doubleScalar a1 = qi1[k];
            qi[k] = c*a0 + s*a1;
            qi1[k] = -s*a0 + c*a1;
        }
    }

    // Row n-1 of R is zero, column n-1 of Q is dropped (free slot).
    // Shift D; the free buffer ends in slot n-1.
    for (label j = 0; j < n - 1; ++j)
    {
        D_[j].swap(D_[j + 1]);
    }

    nHist_ = n - 1;
}


bool Foam::anderson::orthogonalise(const label c)
{
    List<doubleScalar>& v = Q_[c];

    for (label i = 0; i <= c; ++i)
    {
        R(i, c) = 0;
    }

    // Two classical Gram-Schmidt passes (CGS2), one reduction each; the
    // first pass also returns the squared norm before orthogonalisation
    doubleScalar norm0Sqr = 0;

    for (label pass = 0; pass < 2 && c > 0; ++pass)
    {
        List<doubleScalar> h(c + 1, Zero);
        forAll(v, k)
        {
            const doubleScalar vk = v[k];
            for (label j = 0; j < c; ++j)
            {
                h[j] += Q_[j][k]*vk;
            }
            h[c] += vk*vk;
        }
        sumAll(h);

        if (pass == 0)
        {
            norm0Sqr = h[c];
        }

        for (label j = 0; j < c; ++j)
        {
            R(j, c) += h[j];
        }
        forAll(v, k)
        {
            doubleScalar vk = v[k];
            for (label j = 0; j < c; ++j)
            {
                vk -= h[j]*Q_[j][k];
            }
            v[k] = vk;
        }
    }

    List<doubleScalar> nrm(1, Zero);
    forAll(v, k)
    {
        nrm[0] += v[k]*v[k];
    }
    sumAll(nrm);

    if (c == 0)
    {
        norm0Sqr = nrm[0];
    }

    const doubleScalar norm = std::sqrt(nrm[0]);
    const doubleScalar norm0 = std::sqrt(norm0Sqr);

    // Numerical rank: the orthogonal part must keep at least half of the
    // significant digits of the column (sqrt of the double epsilon)
    const doubleScalar rankTol = coupledDefaults::andersonRankTol;

    if (!(norm > rankTol*norm0) || !std::isfinite(norm))
    {
        return false;
    }

    // GUARD: norm > rankTol*norm0 >= 0 checked above
    const doubleScalar rNorm = 1/norm;
    forAll(v, k)
    {
        v[k] *= rNorm;
    }
    R(c, c) = norm;

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::anderson::statusName(const status s)
{
    switch (s)
    {
        case status::applied: return "applied";
        case status::skipped: return "skipped";
        case status::collecting: return "collecting";
        case status::disabled: return "disabled";
    }
    return "unknown";
}


Foam::anderson::status Foam::anderson::apply
(
    volVectorField& U,
    volScalarField& p,
    const vectorField& dU,
    const scalarField& dp,
    const scalar Uref,
    const scalar pref
)
{
    lastMaxAlpha_ = 0;

    if (!enabled_)
    {
        lastStatus_ = status::disabled;
        return lastStatus_;
    }

    allocate();

    // GUARD: reference values floored at VSMALL before inversion
    const doubleScalar sU = 1/static_cast<doubleScalar>(max(Uref, VSMALL));
    const doubleScalar sP = 1/static_cast<doubleScalar>(max(pref, VSMALL));

    // The least-squares norm must be the same for the whole history
    if (hasPrev_ && (sU != sU_ || sP != sP_))
    {
        flush();
    }
    sU_ = sU;
    sP_ = sP;

    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();

    // --- Store the new difference pair, update x_prev, f_prev
    const bool addColumn = hasPrev_;
    if (addColumn && nHist_ == m_)
    {
        removeOldest();
    }
    const label c = nHist_;

    {
        List<doubleScalar>* qPtr = (addColumn ? &Q_[c] : nullptr);
        List<doubleScalar>* dPtr = (addColumn ? &D_[c] : nullptr);

        auto store = [&](const label k, const doubleScalar y,
                         const doubleScalar f, const doubleScalar s)
        {
            // x_k = y_k - f_k (state at which the step was computed)
            const doubleScalar x = y - f;
            if (addColumn)
            {
                const doubleScalar df = f - fPrev_[k];
                (*qPtr)[k] = s*df;
                (*dPtr)[k] = (x - xPrev_[k]) + beta_*df;
            }
            xPrev_[k] = x;
            fPrev_[k] = f;
        };

        forAll(Ui, celli)
        {
            const label k0 = nCmpt*celli;
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                store
                (
                    k0 + d,
                    static_cast<doubleScalar>(Ui[celli][d]),
                    static_cast<doubleScalar>(dU[celli][d]),
                    sU
                );
            }
            store
            (
                k0 + pCmpt,
                static_cast<doubleScalar>(pi[celli]),
                static_cast<doubleScalar>(dp[celli]),
                sP
            );
        }
    }
    hasPrev_ = true;

    if (addColumn && orthogonalise(c))
    {
        ++nHist_;
    }

    if (nHist_ == 0)
    {
        lastStatus_ = status::collecting;
        return lastStatus_;
    }

    const label n = nHist_;

    // --- Least squares: R gamma = Q^T S f_k
    List<doubleScalar> g(n, Zero);
    forAll(Ui, celli)
    {
        const label k0 = nCmpt*celli;
        for (label cmpt = 0; cmpt < nCmpt; ++cmpt)
        {
            const label k = k0 + cmpt;
            const doubleScalar fs = (cmpt == pCmpt ? sP : sU)*fPrev_[k];
            for (label j = 0; j < n; ++j)
            {
                g[j] += Q_[j][k]*fs;
            }
        }
    }
    sumAll(g);

    List<doubleScalar> gamma(n, Zero);
    for (label j = n - 1; j >= 0; --j)
    {
        doubleScalar s = g[j];
        for (label l = j + 1; l < n; ++l)
        {
            s -= R(j, l)*gamma[l];
        }
        // GUARD: R(j, j) > 0 for every accepted column (orthogonalise
        // rank test; Givens downdate keeps the diagonal positive)
        gamma[j] = s/R(j, j);
    }

    doubleScalar maxA = 0;
    for (const doubleScalar a : gamma)
    {
        maxA = std::fmax(maxA, std::fabs(a));
    }
    lastMaxAlpha_ = maxA;
    lastGamma_ = gamma;

    if (!std::isfinite(maxA))
    {
        // Corrupt history: skip and start over
        flush();
        ++nSkipped_;
        lastStatus_ = status::skipped;
        return lastStatus_;
    }

    if (maxA > maxAlpha_)
    {
        // Safeguard (B5): skip, keep the history
        ++nSkipped_;
        lastStatus_ = status::skipped;
        return lastStatus_;
    }

    // --- x_(k+1) = x_k + beta*f_k - sum_j gamma_j D_j
    vectorField& UiRef = U.primitiveFieldRef();
    scalarField& piRef = p.primitiveFieldRef();

    auto extrapolated = [&](const label k)
    {
        doubleScalar z = xPrev_[k] + beta_*fPrev_[k];
        for (label j = 0; j < n; ++j)
        {
            z -= gamma[j]*D_[j][k];
        }
        return z;
    };

    forAll(UiRef, celli)
    {
        const label k0 = nCmpt*celli;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            UiRef[celli][d] = static_cast<scalar>(extrapolated(k0 + d));
        }
        piRef[celli] = static_cast<scalar>(extrapolated(k0 + pCmpt));
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    ++nApplied_;
    lastStatus_ = status::applied;
    return lastStatus_;
}


void Foam::anderson::reject(volVectorField& U, volScalarField& p)
{
    if (!hasPrev_)
    {
        FatalErrorInFunction
            << "anderson: no un-extrapolated state to restore"
            << exit(FatalError);
    }

    // y_k = x_k + f_k
    vectorField& Ui = U.primitiveFieldRef();
    scalarField& pi = p.primitiveFieldRef();
    forAll(Ui, celli)
    {
        const label k0 = nCmpt*celli;
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            Ui[celli][d] =
                static_cast<scalar>(xPrev_[k0 + d] + fPrev_[k0 + d]);
        }
        pi[celli] =
            static_cast<scalar>(xPrev_[k0 + pCmpt] + fPrev_[k0 + pCmpt]);
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    ++nRejected_;
    flush();
}


void Foam::anderson::flush(const char* reason)
{
    if (hasPrev_ || nHist_ > 0)
    {
        ++nFlushed_;
        lastFlushReason_ = reason;
    }
    nHist_ = 0;
    hasPrev_ = false;
}


Foam::doubleScalar Foam::anderson::conditionEstimate() const
{
    if (nHist_ == 0)
    {
        return 0;
    }
    doubleScalar mx = 0;
    doubleScalar mn = GREAT;
    for (label j = 0; j < nHist_; ++j)
    {
        const doubleScalar r = std::fabs(R(j, j));
        mx = std::fmax(mx, r);
        mn = std::fmin(mn, r);
    }
    // GUARD: accepted columns have R(j, j) > 0 (rank test)
    return mx/std::fmax(mn, doubleScalarVSMALL);
}


void Foam::anderson::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("enabled", enabled_);
    d.add("m", m_);
    d.add("beta", beta_);
    d.add("maxAlpha", maxAlpha_);
    dict.add("anderson", d);
}


// ************************************************************************* //
