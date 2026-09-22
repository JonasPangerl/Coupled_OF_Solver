/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockSimplePrecon.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(blockSimplePrecon, 0);
    addToRunTimeSelectionTable
    (
        blockPreconditioner,
        blockSimplePrecon,
        dictionary
    );
}


namespace
{
    // Index of the pressure unknown in a block
    constexpr Foam::label pIdx = Foam::blockDim - 1;

    // Copy a block keeping u-u and scaled p-p entries (cross blocks zero)
    inline void decouple
    (
        const Foam::blockScalar* src,
        Foam::blockScalar* dst,
        const Foam::blockScalar sScale
    )
    {
        using namespace Foam;
        for (label r = 0; r < blockDim; ++r)
        {
            for (label c = 0; c < blockDim; ++c)
            {
                const label i = r*blockDim + c;
                if (r < pIdx && c < pIdx)
                {
                    dst[i] = src[i];
                }
                else if (r == pIdx && c == pIdx)
                {
                    dst[i] = sScale*src[i];
                }
                else
                {
                    dst[i] = 0;
                }
            }
        }
    }
}


Foam::blockSimplePrecon::blockSimplePrecon
(
    const blockSolver& solver,
    const dictionary& dict
)
:
    blockPreconditioner(solver),
    Md_(solver.matrix().mesh()),
    gamg_(Md_, dict.subOrEmptyDict("blockGAMG")),
    schurScale_
    (
        dict.getOrDefault<doubleScalar>
        (
            "schurScale",
            coupledDefaults::schurScale
        )
    ),
    pivotGuard_
    (
        dict.getOrDefault<doubleScalar>
        (
            "pivotGuard",
            coupledDefaults::pivotGuard
        )
    ),
    sequential_(false),
    rL_(),
    t_(),
    v_(),
    y_(),
    y2_()
{
    const word mode(dict.getOrDefault<word>
    (
        "simpleMode",
        word(coupledDefaults::simpleMode)
    ));
    if (mode == "sequential")
    {
        sequential_ = true;
    }
    else if (mode != "single")
    {
        FatalIOErrorInFunction(dict)
            << "simpleMode " << mode << ": valid single sequential"
            << exit(FatalIOError);
    }
    if (schurScale_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "schurScale " << schurScale_ << " must be > 0"
            << exit(FatalIOError);
    }
    const word solverType(dict.getOrDefault<word>("solver", word::null));
    if (solverType != "blockFGMRES")
    {
        FatalIOErrorInFunction(dict)
            << "blockSimple is a variable preconditioner: use blockFGMRES"
            << exit(FatalIOError);
    }
    gamg_.writeStats(Info);
}


void Foam::blockSimplePrecon::update()
{
    const blockLduMatrix4& A = matrix();
    const label nCells = A.nCells();
    const blockScalar s = narrow(schurScale_);

    for (label c = 0; c < nCells; ++c)
    {
        decouple(A.diagBlock(c), Md_.diagBlock(c), s);
    }
    const label nF = A.nFaces();
    for (label f = 0; f < nF; ++f)
    {
        decouple
        (
            A.upper().cdata() + f*blockSize,
            Md_.upper().data() + f*blockSize,
            s
        );
        decouple
        (
            A.lower().cdata() + f*blockSize,
            Md_.lower().data() + f*blockSize,
            s
        );
    }
    forAll(A.interfaces(), i)
    {
        const blockScalarList& src = A.interfaceCoeffs(i);
        blockScalarList& dst = Md_.interfaceCoeffs(i);
        const label nb = src.size()/blockSize;
        for (label b = 0; b < nb; ++b)
        {
            decouple(src.cdata() + b*blockSize, dst.data() + b*blockSize, s);
        }
    }
    Md_.markUpdated();
    gamg_.update();

    // 1/diag of the u-u block
    rL_.resize_nocopy(blockDim*nCells);
    const reduceScalar guard = pivotGuard_;
    for (label c = 0; c < nCells; ++c)
    {
        for (label k = 0; k < pIdx; ++k)
        {
            const reduceScalar d = toDouble(A.diagBlock(c)[k*blockDim + k]);
            // GUARD: pivot guard as block4Ops
            const reduceScalar dg =
                (std::abs(d) < guard ? (d < 0 ? -guard : guard) : d);
            rL_[c*blockDim + k] = narrow(1.0/dg);
        }
        rL_[c*blockDim + pIdx] = 0;
    }

    const label n = A.nRows();
    t_.resize_nocopy(n);
    v_.resize_nocopy(n);
    y_.resize_nocopy(n);
    y2_.resize_nocopy(n);
}


void Foam::blockSimplePrecon::precondition
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    const blockLduMatrix4& A = matrix();
    const label nCells = A.nCells();

    if (sequential_)
    {
        // y_u = A^-1 r_u (one cycle on the momentum part)
        t_ = r;
        for (label c = 0; c < nCells; ++c)
        {
            t_[c*blockDim + pIdx] = 0;
        }
        gamg_.apply(y_, t_);
        for (label c = 0; c < nCells; ++c)
        {
            y_[c*blockDim + pIdx] = 0;
        }
    }
    else
    {
        // Predictor y_u = L^-1 r_u
        for (label c = 0; c < nCells; ++c)
        {
            for (label k = 0; k < pIdx; ++k)
            {
                const label i = c*blockDim + k;
                y_[i] = rL_[i]*r[i];
            }
            y_[c*blockDim + pIdx] = 0;
        }
    }

    // r_p - B y_u  (p rows of A (y_u, 0))
    A.Amul(v_, y_);
    t_ = r;
    for (label c = 0; c < nCells; ++c)
    {
        t_[c*blockDim + pIdx] -= v_[c*blockDim + pIdx];
        if (sequential_)
        {
            for (label k = 0; k < pIdx; ++k)
            {
                t_[c*blockDim + k] = 0;
            }
        }
    }
    gamg_.apply(y2_, t_);

    if (sequential_)
    {
        // y2 carries the p part only; the u part is from the first cycle
        for (label c = 0; c < nCells; ++c)
        {
            for (label k = 0; k < pIdx; ++k)
            {
                y2_[c*blockDim + k] = y_[c*blockDim + k];
            }
        }
    }

    // u = y_u - L^-1 G y_p  (u rows of A (0, y_p))
    for (label c = 0; c < nCells; ++c)
    {
        for (label k = 0; k < pIdx; ++k)
        {
            y_[c*blockDim + k] = 0;
        }
        y_[c*blockDim + pIdx] = y2_[c*blockDim + pIdx];
    }
    A.Amul(v_, y_);
    for (label c = 0; c < nCells; ++c)
    {
        for (label k = 0; k < pIdx; ++k)
        {
            const label i = c*blockDim + k;
            w[i] = y2_[i] - rL_[i]*v_[i];
        }
        w[c*blockDim + pIdx] = y2_[c*blockDim + pIdx];
    }
}


void Foam::blockSimplePrecon::writeSettings(dictionary& dict) const
{
    dict.add("type", type());
    dict.add("schurScale", schurScale_);
    dict.add("simpleMode", word(sequential_ ? "sequential" : "single"));
    dictionary gd;
    gamg_.writeSettings(gd);
    dict.add("blockGAMG", gd);
}


// ************************************************************************* //
