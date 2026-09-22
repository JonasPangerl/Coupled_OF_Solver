/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-kernelBandwidth

Description
    Phase A gate of spec 6.5.6 (amendment B11): achieved memory bandwidth of
    the hot-loop kernels against the STREAM triad of the same machine, one
    thread, same allocation and alignment (alignedList, 64 B).

    System: synthetic structured nx^3 grid (default 171^3 = 5.0 M cells),
    7-point connectivity, LDU face order (faces sorted by owner, then by
    neighbour), dense 4x4 float blocks on the diagonal and on every face
    (blockLduMatrix4 on an lduPrimitiveMesh, no interfaces).

    Kernels and bytes counted (every operand read once, written once; the
    write-allocate read of a written array is not counted, as in STREAM):
        triad         a = b + s c              3 x 4 B per entry
        Amul          Ax = A x                 diag 64 B + x 16 B + Ax 16 B
                                               per cell, upper + lower 128 B
                                               + addressing 2 x 4 B per face
        axpy_dot      y += a x; <y,y>, <y,z>   4 x 4 B per entry
        update_residual_norm, mgs_axpy_dot and dot are reported for
        information (not gated).

    Every kernel runs `reps` times after one warm-up call; the best time is
    used (STREAM convention), the median is reported too.

    Pass: Amul >= 60 % and axpy_dot >= 80 % of the triad bandwidth
    (coupledDefaults::gateAmulFraction, gateAxpyDotFraction).

Usage
    Test-kernelBandwidth [-nx N] [-reps R] [-json <file>]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "blockLduMatrix4.H"
#include "blockKernels.H"
#include "alignedList.H"
#include "coupledDefaults.H"
#include "jsonWriter.H"
#include "lduPrimitiveMesh.H"
#include <chrono>
#include <algorithm>
#include <cstdint>
#include <vector>

using namespace Foam;

// Test parameters
static constexpr label defaultNx = 171;       // 171^3 = 5.0 M cells (6.5.6)
static constexpr label defaultReps = 10;
static constexpr label nWarmup = 1;
static constexpr doubleScalar bytesPerGB = 1e9;
static constexpr blockScalar triadScalar = 3;
// Small axpy coefficient: y stays bounded over the repetitions (the sign
// alternates) and all values stay normal floats
static constexpr blockScalar axpyCoeff = 1e-3f;
// Diagonal dominance of the synthetic blocks: |off-diagonal entries| <= 1
static constexpr blockScalar diagShift = 32;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- Deterministic values in [0.5, 1) with a pseudo-random sign (LCG)
class lcg
{
    std::uint64_t s_;

public:

    explicit lcg(const std::uint64_t seed) : s_(seed) {}

    blockScalar next()
    {
        // Knuth MMIX constants
        s_ = s_*6364136223846793005ULL + 1442695040888963407ULL;
        const std::uint32_t top = static_cast<std::uint32_t>(s_ >> 40);
        const blockScalar u =
            static_cast<blockScalar>(top & 0xFFFFu)/65536.0f;
        const blockScalar v = 0.5f + 0.5f*u;
        return (top & 0x10000u) ? -v : v;
    }
};


struct timing
{
    doubleScalar best;
    doubleScalar median;
};


template<class Fn>
timing timeKernel(const label reps, Fn&& fn)
{
    for (label i = 0; i < nWarmup; ++i)
    {
        fn(i);
    }
    std::vector<doubleScalar> t;
    t.reserve(static_cast<std::size_t>(reps));
    for (label i = 0; i < reps; ++i)
    {
        const auto t0 = std::chrono::steady_clock::now();
        fn(i);
        const auto t1 = std::chrono::steady_clock::now();
        t.push_back(std::chrono::duration<doubleScalar>(t1 - t0).count());
    }
    std::sort(t.begin(), t.end());
    return {t.front(), t[t.size()/2]};
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("nx", "N", "Cells per direction (171 -> 5.0 M)");
    argList::addOption("reps", "R", "Timed repetitions per kernel (10)");
    argList::addOption("json", "file", "Write results as JSON");
    argList args(argc, argv);

    const label nx = args.getOrDefault<label>("nx", defaultNx);
    const label reps = max(label(1), args.getOrDefault<label>("reps", defaultReps));
    const label nCells = nx*nx*nx;

    // Structured 7-point connectivity in LDU order
    labelList lower;
    labelList upper;
    {
        DynamicList<label> l(3*nCells);
        DynamicList<label> u(3*nCells);
        for (label k = 0; k < nx; ++k)
        {
            for (label j = 0; j < nx; ++j)
            {
                for (label i = 0; i < nx; ++i)
                {
                    const label c = i + nx*(j + nx*k);
                    if (i + 1 < nx) { l.push_back(c); u.push_back(c + 1); }
                    if (j + 1 < nx) { l.push_back(c); u.push_back(c + nx); }
                    if (k + 1 < nx) { l.push_back(c); u.push_back(c + nx*nx); }
                }
            }
        }
        lower.transfer(l);
        upper.transfer(u);
    }
    const label nFaces = lower.size();

    lduPrimitiveMesh mesh(nCells, lower, upper, UPstream::worldComm, true);
    blockLduMatrix4 A(mesh);

    {
        lcg rnd(1);
        for (blockScalar& v : A.upper()) v = rnd.next();
        for (blockScalar& v : A.lower()) v = rnd.next();
        blockScalar* d = A.diag().data();
        for (label c = 0; c < nCells; ++c)
        {
            for (label e = 0; e < blockSize; ++e)
            {
                d[c*blockSize + e] = rnd.next();
            }
            for (label r = 0; r < blockDim; ++r)
            {
                d[c*blockSize + r*blockDim + r] += diagShift;
            }
        }
    }

    const label n = blockDim*nCells;

    alignedList<blockScalar> x(n);
    alignedList<blockScalar> Ax(n);
    alignedList<blockScalar> y(n);
    alignedList<blockScalar> z(n);
    {
        lcg rnd(2);
        for (label i = 0; i < n; ++i)
        {
            x[i] = rnd.next();
            y[i] = rnd.next();
            z[i] = rnd.next();
        }
        Ax.setZero();
    }

    Info<< "Test-kernelBandwidth: nx " << nx << ", cells " << nCells
        << ", faces " << nFaces << ", unknowns " << n
        << ", reps " << reps << nl
        << "  matrix storage " << doubleScalar(A.storageBytes())/bytesPerGB
        << " GB" << endl;

    const bool alignedMatrix =
        alignedList<blockScalar>::isAligned(A.diag().cdata())
     && alignedList<blockScalar>::isAligned(A.upper().cdata())
     && alignedList<blockScalar>::isAligned(A.lower().cdata());
    const bool alignedVectors =
        alignedList<blockScalar>::isAligned(x.cdata())
     && alignedList<blockScalar>::isAligned(y.cdata());

    // Keeps the reductions alive
    volatile reduceScalar sink = 0;

    // STREAM triad a = b + s c on three arrays of n entries (the triad
    // reuses Ax, x, z: same allocation and alignment as the kernels)
    const timing tTriad = timeKernel
    (
        reps,
        [&](label)
        {
            blockKernels::triad(n, triadScalar, x.cdata(), z.cdata(), Ax.data());
        }
    );
    sink = sink + toDouble(Ax[n/2]);

    // Amul
    const timing tAmul = timeKernel
    (
        reps,
        [&](label) { A.Amul(Ax, x); }
    );
    sink = sink + toDouble(Ax[n/2]);

    // axpy_dot (alternating sign keeps y bounded)
    const timing tAxpyDot = timeKernel
    (
        reps,
        [&](label rep)
        {
            const blockScalar a = (rep % 2 ? -axpyCoeff : axpyCoeff);
            const blockKernels::sumPair s =
                blockKernels::axpy_dot(n, a, x.cdata(), y.data(), z.cdata());
            sink = sink + s.first + s.second;
        }
    );

    // update_residual_norm
    const timing tUpdRes = timeKernel
    (
        reps,
        [&](label rep)
        {
            const blockScalar a = (rep % 2 ? -axpyCoeff : axpyCoeff);
            sink = sink
              + blockKernels::update_residual_norm(n, a, x.cdata(), y.data());
        }
    );

    // mgs_axpy_dot
    const timing tMgs = timeKernel
    (
        reps,
        [&](label rep)
        {
            const blockScalar a = (rep % 2 ? -axpyCoeff : axpyCoeff);
            sink = sink
              + blockKernels::mgs_axpy_dot
                (
                    n, a, x.cdata(), y.data(), z.cdata()
                );
        }
    );

    // dot
    const timing tDot = timeKernel
    (
        reps,
        [&](label)
        {
            sink = sink + blockKernels::dot(n, x.cdata(), z.cdata());
        }
    );

    const doubleScalar F = sizeof(blockScalar);
    const doubleScalar L = sizeof(label);
    const doubleScalar nD = doubleScalar(n);

    const doubleScalar bTriad = 3*F*nD;
    const doubleScalar bAmul =
        doubleScalar(nCells)*(blockSize + 2*blockDim)*F
      + doubleScalar(nFaces)*(2*blockSize*F + 2*L);
    const doubleScalar bAxpyDot = 4*F*nD;
    const doubleScalar bUpdRes = 3*F*nD;
    const doubleScalar bMgs = 4*F*nD;
    const doubleScalar bDot = 2*F*nD;

    auto gbs = [](const doubleScalar bytes, const doubleScalar t)
    {
        return bytes/t/bytesPerGB;  // GUARD: t > 0 (measured wall time)
    };

    const doubleScalar gTriad = gbs(bTriad, tTriad.best);
    const doubleScalar gAmul = gbs(bAmul, tAmul.best);
    const doubleScalar gAxpyDot = gbs(bAxpyDot, tAxpyDot.best);
    const doubleScalar gUpdRes = gbs(bUpdRes, tUpdRes.best);
    const doubleScalar gMgs = gbs(bMgs, tMgs.best);
    const doubleScalar gDot = gbs(bDot, tDot.best);

    const doubleScalar fAmul = gAmul/gTriad;
    const doubleScalar fAxpyDot = gAxpyDot/gTriad;

    const bool passAmul = fAmul >= coupledDefaults::gateAmulFraction;
    const bool passAxpyDot = fAxpyDot >= coupledDefaults::gateAxpyDotFraction;
    const bool pass = passAmul && passAxpyDot;

    auto line = [&](const char* name, const doubleScalar g, const timing& t)
    {
        Info<< "  " << name << ": " << g << " GB/s ("
            << 100*g/gTriad << " % of triad), best " << t.best
            << " s, median " << t.median << " s" << endl;
    };
    line("triad               ", gTriad, tTriad);
    line("Amul                ", gAmul, tAmul);
    line("axpy_dot            ", gAxpyDot, tAxpyDot);
    line("update_residual_norm", gUpdRes, tUpdRes);
    line("mgs_axpy_dot        ", gMgs, tMgs);
    line("dot                 ", gDot, tDot);
    Info<< "  alignment 64 B: matrix " << alignedMatrix
        << ", vectors " << alignedVectors << nl
        << "  Amul gate (>= " << 100*coupledDefaults::gateAmulFraction
        << " %): " << (passAmul ? "PASS" : "FAIL") << nl
        << "  axpy_dot gate (>= " << 100*coupledDefaults::gateAxpyDotFraction
        << " %): " << (passAxpyDot ? "PASS" : "FAIL") << nl
        << "  checksum " << reduceScalar(sink) << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-kernelBandwidth");
        j.add("gate", "phase_A (spec 6.5.6)");
        j.add("nx", nx);
        j.add("nCells", nCells);
        j.add("nFaces", nFaces);
        j.add("nUnknowns", n);
        j.add("reps", reps);
        j.add("threads", label(1));
        j.add("alignedMatrix64", alignedMatrix);
        j.add("alignedVectors64", alignedVectors);
        j.add("triad_GBps", gTriad);
        j.add("triad_bytes", bTriad);
        j.add("triad_tBest", tTriad.best);
        j.add("Amul_GBps", gAmul);
        j.add("Amul_bytes", bAmul);
        j.add("Amul_tBest", tAmul.best);
        j.add("Amul_tMedian", tAmul.median);
        j.add("Amul_fractionOfTriad", fAmul);
        j.add("axpy_dot_GBps", gAxpyDot);
        j.add("axpy_dot_bytes", bAxpyDot);
        j.add("axpy_dot_tBest", tAxpyDot.best);
        j.add("axpy_dot_tMedian", tAxpyDot.median);
        j.add("axpy_dot_fractionOfTriad", fAxpyDot);
        j.add("update_residual_norm_GBps", gUpdRes);
        j.add("mgs_axpy_dot_GBps", gMgs);
        j.add("dot_GBps", gDot);
        j.add("gateAmulFraction", coupledDefaults::gateAmulFraction);
        j.add("gateAxpyDotFraction", coupledDefaults::gateAxpyDotFraction);
        j.add("passAmul", passAmul);
        j.add("passAxpyDot", passAxpyDot);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
