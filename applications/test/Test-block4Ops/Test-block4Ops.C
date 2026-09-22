/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-block4Ops

Description
    Unit test of the dense 4x4 block kernels (spec 6.4):

    1. invert 1000 random well-conditioned 4x4 blocks (float storage),
       ||A A^-1 - I||_inf < 1e-5 for every block;
    2. a singular block (zero, and rank-deficient) triggers the pivot guard
       and produces no NaN/Inf;
    3. matVec / matMul consistency against double reference products.
    4. 3x3 kernels of the tensorial Rhie-Chow D (amendment C2), double:
       invert3 on random well-conditioned (non-symmetric) matrices,
       ||A A^-1 - I||_inf < 1e-12; pseudoInverse3 of the same matrices
       equals the inverse to 1e-10; a rank-2 matrix and the zero matrix are
       rejected by the invert3 determinant guard, and their pseudo-inverse
       satisfies the Penrose conditions A A+ A = A, A+ A A+ = A+ and
       symmetric A A+, A+ A to 1e-10 (relative), with the expected rank.

    Well-conditioned: A = R + s I with R uniform in [-1, 1] and s = 4, which
    makes A strictly diagonally dominant (condition number O(10)).

Usage
    Test-block4Ops [-json <file>] [-nBlocks <n>] [-seed <n>]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Random.H"
#include "block4Ops.H"
#include "coupledDefaults.H"
#include "jsonWriter.H"
#include <cmath>

using namespace Foam;

// Test parameters (not solver parameters; fixed by spec 6.4)
static constexpr label defaultNBlocks = 1000;
static constexpr doubleScalar tolInverse3 = 1e-12;
static constexpr doubleScalar tolPinv3 = 1e-10;

// 3x3 helpers of test 4 (row-major, double)
static void mul3(const doubleScalar* A, const doubleScalar* B, doubleScalar* C)
{
    for (label r = 0; r < 3; ++r)
    {
        for (label c = 0; c < 3; ++c)
        {
            doubleScalar s = 0;
            for (label k = 0; k < 3; ++k)
            {
                s += A[r*3 + k]*B[k*3 + c];
            }
            C[r*3 + c] = s;
        }
    }
}

static doubleScalar maxAbs3(const doubleScalar* A)
{
    doubleScalar m = 0;
    for (label i = 0; i < 9; ++i) m = std::max(m, std::abs(A[i]));
    return m;
}

//- Largest relative violation of the four Penrose conditions
static doubleScalar penroseError(const doubleScalar* A, const doubleScalar* P)
{
    doubleScalar AP[9], PA[9], APA[9], PAP[9];
    mul3(A, P, AP);
    mul3(P, A, PA);
    mul3(AP, A, APA);
    mul3(PA, P, PAP);
    const doubleScalar nA = std::max(maxAbs3(A), doubleScalar(1e-300));
    const doubleScalar nP = std::max(maxAbs3(P), doubleScalar(1e-300));
    doubleScalar e = 0;
    for (label r = 0; r < 3; ++r)
    {
        for (label c = 0; c < 3; ++c)
        {
            const label i = r*3 + c;
            const label t = c*3 + r;
            e = std::max(e, std::abs(APA[i] - A[i])/nA);
            e = std::max(e, std::abs(PAP[i] - P[i])/nP);
            e = std::max(e, std::abs(AP[i] - AP[t]));
            e = std::max(e, std::abs(PA[i] - PA[t]));
        }
    }
    return e;
}
static constexpr doubleScalar tolInverse = 1e-5;
static constexpr doubleScalar diagShift = 4;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static doubleScalar invError(const blockScalar* A, const blockScalar* Ai)
{
    doubleScalar err = 0;
    for (label r = 0; r < blockDim; ++r)
    {
        for (label c = 0; c < blockDim; ++c)
        {
            doubleScalar s = 0;
            for (label k = 0; k < blockDim; ++k)
            {
                s += toDouble(A[r*blockDim + k])*toDouble(Ai[k*blockDim + c]);
            }
            s -= (r == c ? 1.0 : 0.0);
            err = std::max(err, std::abs(s));
        }
    }
    return err;
}


static bool allFinite(const blockScalar* A)
{
    for (label i = 0; i < blockSize; ++i)
    {
        if (!std::isfinite(A[i])) return false;
    }
    return true;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("json", "file", "Write results as JSON");
    argList::addOption("nBlocks", "n", "Number of random blocks (1000)");
    argList::addOption("seed", "n", "Random seed (1)");
    argList args(argc, argv);

    const label nBlocks = args.getOrDefault<label>("nBlocks", defaultNBlocks);
    Random rnd(args.getOrDefault<label>("seed", 1));

    const doubleScalar pivotGuard = coupledDefaults::pivotGuard;

    // 1. Random well-conditioned blocks
    doubleScalar maxErr = 0;
    label nFail = 0;
    label nGuardedRandom = 0;
    blockScalar A[blockSize];
    blockScalar Ai[blockSize];

    for (label b = 0; b < nBlocks; ++b)
    {
        for (label i = 0; i < blockSize; ++i)
        {
            A[i] = narrow(2*doubleScalar(rnd.sample01<scalar>()) - 1);
        }
        for (label i = 0; i < blockDim; ++i)
        {
            A[i*blockDim + i] += narrow(diagShift);
        }

        nGuardedRandom += block4Ops::invert(A, Ai, pivotGuard);
        const doubleScalar e = invError(A, Ai);
        maxErr = std::max(maxErr, e);
        if (!(e < tolInverse)) ++nFail;
    }

    Info<< "random blocks: " << nBlocks
        << "  max ||A A^-1 - I||_inf = " << maxErr
        << "  failures " << nFail
        << "  guarded pivots " << nGuardedRandom << endl;

    // 2a. Zero block
    block4Ops::zero(A);
    const label nGuardZero = block4Ops::invert(A, Ai, pivotGuard);
    const bool finiteZero = allFinite(Ai);

    // 2b. Exactly singular block with a regular momentum part: random
    //     well-conditioned 3x3 u-v-w block, continuity row and pressure
    //     column zero (a cell without p-U coupling). A rank deficiency built
    //     from float sums would leave a pivot of O(1e-8) instead of 0.
    for (label i = 0; i < blockSize; ++i)
    {
        A[i] = narrow(2*doubleScalar(rnd.sample01<scalar>()) - 1);
    }
    for (label i = 0; i < blockP; ++i)
    {
        A[i*blockDim + i] += narrow(diagShift);
    }
    for (label c = 0; c < blockDim; ++c)
    {
        A[blockP*blockDim + c] = 0;
        A[c*blockDim + blockP] = 0;
    }
    const label nGuardRank = block4Ops::invert(A, Ai, pivotGuard);
    const bool finiteRank = allFinite(Ai);

    Info<< "zero block: guarded pivots " << nGuardZero
        << " finite " << finiteZero << nl
        << "singular block: guarded pivots " << nGuardRank
        << " finite " << finiteRank << endl;

    // 3. matVec against a double reference
    doubleScalar maxMatVecErr = 0;
    for (label b = 0; b < nBlocks; ++b)
    {
        blockScalar x[blockDim];
        blockScalar y[blockDim];
        for (label i = 0; i < blockSize; ++i)
        {
            A[i] = narrow(2*doubleScalar(rnd.sample01<scalar>()) - 1);
        }
        for (label i = 0; i < blockDim; ++i)
        {
            x[i] = narrow(2*doubleScalar(rnd.sample01<scalar>()) - 1);
        }
        block4Ops::matVec(A, x, y);
        for (label r = 0; r < blockDim; ++r)
        {
            doubleScalar ref = 0;
            for (label c = 0; c < blockDim; ++c)
            {
                ref += toDouble(A[r*blockDim + c])*toDouble(x[c]);
            }
            maxMatVecErr = std::max(maxMatVecErr, std::abs(ref - toDouble(y[r])));
        }
    }
    Info<< "matVec max abs error vs double reference = " << maxMatVecErr
        << endl;

    // 4. 3x3 kernels (C2)
    const doubleScalar detTol = coupledDefaults::rhieChowDetRelTol;
    const doubleScalar pinvTol = coupledDefaults::rhieChowPinvRelTol;
    const label sweeps = coupledDefaults::rhieChowPinvMaxSweeps;
    doubleScalar maxErr3 = 0;
    doubleScalar maxPinvVsInv = 0;
    label nRejected3 = 0;
    for (label b = 0; b < nBlocks; ++b)
    {
        doubleScalar A3[9], I3[9], P3[9], AI[9];
        for (label i = 0; i < 9; ++i)
        {
            A3[i] = 2*rnd.sample01<doubleScalar>() - 1;
        }
        for (label i = 0; i < 3; ++i)
        {
            A3[i*3 + i] += diagShift;
        }
        if (!block4Ops::invert3(A3, I3, detTol))
        {
            ++nRejected3;
            continue;
        }
        mul3(A3, I3, AI);
        for (label r = 0; r < 3; ++r)
        {
            for (label c = 0; c < 3; ++c)
            {
                maxErr3 = std::max
                (
                    maxErr3, std::abs(AI[r*3 + c] - (r == c ? 1.0 : 0.0))
                );
            }
        }
        block4Ops::pseudoInverse3(A3, P3, pinvTol, sweeps);
        for (label i = 0; i < 9; ++i)
        {
            maxPinvVsInv = std::max
            (
                maxPinvVsInv, std::abs(P3[i] - I3[i])/maxAbs3(I3)
            );
        }
    }

    // Rank-2 matrix: third row = row 0 + 2 row 1 (exact in double for
    // these dyadic values), non-symmetric
    const doubleScalar R2[9] =
    {
        2.0, -1.0, 0.5,
        0.25, 3.0, -1.5,
        2.5, 5.0, -2.5
    };
    doubleScalar tmp3[9], P2[9];
    const bool rejectedRank2 = !block4Ops::invert3(R2, tmp3, detTol);
    const label rank2 = block4Ops::pseudoInverse3(R2, P2, pinvTol, sweeps);
    const doubleScalar penrose2 = penroseError(R2, P2);

    const doubleScalar Z3[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    doubleScalar PZ[9];
    const bool rejectedZero3 = !block4Ops::invert3(Z3, tmp3, detTol);
    const label rankZero = block4Ops::pseudoInverse3(Z3, PZ, pinvTol, sweeps);
    const bool zeroPinvZero = (maxAbs3(PZ) == 0);

    Info<< "3x3 invert3: max ||A A^-1 - I||_inf = " << maxErr3
        << "  rejected " << nRejected3
        << "  max |pinv - inv|/|inv| = " << maxPinvVsInv << nl
        << "3x3 rank-2: rejected " << rejectedRank2 << " rank " << rank2
        << " Penrose error " << penrose2 << nl
        << "3x3 zero: rejected " << rejectedZero3 << " rank " << rankZero
        << " pinv zero " << zeroPinvZero << endl;

    const bool pass3 =
        nRejected3 == 0
     && maxErr3 < tolInverse3
     && maxPinvVsInv < tolPinv3
     && rejectedRank2 && rank2 == 2 && penrose2 < tolPinv3
     && rejectedZero3 && rankZero == 0 && zeroPinvZero;

    const bool pass =
        nFail == 0
     && nGuardZero > 0 && finiteZero
     && nGuardRank > 0 && finiteRank
     && maxMatVecErr < tolInverse
     && pass3;

    Info<< (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-block4Ops");
        j.add("nBlocks", nBlocks);
        j.add("maxInverseError", maxErr);
        j.add("tolerance", tolInverse);
        j.add("nFail", nFail);
        j.add("nGuardedRandom", nGuardedRandom);
        j.add("zeroBlockGuarded", nGuardZero);
        j.add("zeroBlockFinite", finiteZero);
        j.add("rankDeficientGuarded", nGuardRank);
        j.add("rankDeficientFinite", finiteRank);
        j.add("maxMatVecError", maxMatVecErr);
        j.add("maxInverse3Error", maxErr3);
        j.add("maxPinv3VsInverse", maxPinvVsInv);
        j.add("rank2Rejected", rejectedRank2);
        j.add("rank2PinvRank", rank2);
        j.add("rank2PenroseError", penrose2);
        j.add("zero3Rejected", rejectedZero3);
        j.add("zero3PinvRank", rankZero);
        j.add("pass3x3", pass3);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
