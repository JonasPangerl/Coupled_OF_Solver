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
            A[i] = narrow(2*rnd.sample01<doubleScalar>() - 1);
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
        A[i] = narrow(2*rnd.sample01<doubleScalar>() - 1);
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
            A[i] = narrow(2*rnd.sample01<doubleScalar>() - 1);
        }
        for (label i = 0; i < blockDim; ++i)
        {
            x[i] = narrow(2*rnd.sample01<doubleScalar>() - 1);
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

    const bool pass =
        nFail == 0
     && nGuardZero > 0 && finiteZero
     && nGuardRank > 0 && finiteRank
     && maxMatVecErr < tolInverse;

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
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
