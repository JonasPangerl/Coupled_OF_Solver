/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockLduMatrix4.H"
#include "block4Ops.H"
#include "lduMesh.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockLduMatrix4::blockLduMatrix4
(
    const lduAddressing& addr,
    const lduInterfacePtrsList& interfaces,
    const label comm
)
:
    addr_(addr),
    meshPtr_(nullptr),
    version_(0),
    comm_(comm),
    nCells_(addr.size()),
    nFaces_(addr.lowerAddr().size()),
    diag_(blockSize*nCells_, Zero),
    upper_(blockSize*nFaces_, Zero),
    lower_(blockSize*nFaces_, Zero),
    source_(blockDim*nCells_, Zero),
    interfaces_(),
    interfaceCoeffs_(),
    residualAcc_(),
    residualHi_(),
    residualLo_()
{
    label nBlock = 0;
    forAll(interfaces, i)
    {
        if (interfaces.set(i) && blockLduInterface::isBlockCoupled(interfaces[i]))
        {
            ++nBlock;
        }
    }

    interfaces_.resize(nBlock);
    interfaceCoeffs_.resize(nBlock);

    nBlock = 0;
    forAll(interfaces, i)
    {
        if (interfaces.set(i) && blockLduInterface::isBlockCoupled(interfaces[i]))
        {
            interfaces_.set(nBlock, new blockLduInterface(interfaces[i], i));
            interfaceCoeffs_[nBlock].resize
            (
                blockSize*interfaces[i].faceCells().size(),
                Zero
            );
            ++nBlock;
        }
    }
}


Foam::blockLduMatrix4::blockLduMatrix4(const lduMesh& mesh)
:
    blockLduMatrix4(mesh.lduAddr(), mesh.interfaces(), mesh.comm())
{
    meshPtr_ = &mesh;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::blockLduMatrix4::mesh() const
{
    if (!meshPtr_)
    {
        FatalErrorInFunction
            << "Matrix was constructed from addressing only, no lduMesh"
            << abort(FatalError);
    }
    return *meshPtr_;
}


Foam::label Foam::blockLduMatrix4::nnzBlocks() const
{
    label n = nCells_ + 2*nFaces_;
    forAll(interfaces_, i)
    {
        n += interfaces_[i].size();
    }
    return n;
}


std::size_t Foam::blockLduMatrix4::storageBytes() const
{
    std::size_t n =
        diag_.size() + upper_.size() + lower_.size() + source_.size();
    forAll(interfaceCoeffs_, i)
    {
        n += interfaceCoeffs_[i].size();
    }
    return n*sizeof(blockScalar);
}


void Foam::blockLduMatrix4::clear()
{
    diag_ = Zero;
    upper_ = Zero;
    lower_ = Zero;
    source_ = Zero;
    forAll(interfaceCoeffs_, i)
    {
        interfaceCoeffs_[i] = Zero;
    }
}


void Foam::blockLduMatrix4::scaleRow(const label r, const blockScalar s)
{
    auto scaleBlocks = [r, s](blockScalarList& a)
    {
        const label nBlocks = a.size()/blockSize;
        for (label b = 0; b < nBlocks; ++b)
        {
            blockScalar* row = a.data() + b*blockSize + r*blockDim;
            for (label c = 0; c < blockDim; ++c)
            {
                row[c] *= s;
            }
        }
    };

    scaleBlocks(diag_);
    scaleBlocks(upper_);
    scaleBlocks(lower_);
    forAll(interfaceCoeffs_, i)
    {
        scaleBlocks(interfaceCoeffs_[i]);
    }

    for (label celli = 0; celli < nCells_; ++celli)
    {
        source_[celli*blockDim + r] *= s;
    }
}


void Foam::blockLduMatrix4::initInterfaces(const blockScalarUList& x) const
{
    forAll(interfaces_, i)
    {
        interfaces_[i].initExchange(x);
    }
}


void Foam::blockLduMatrix4::updateInterfaces
(
    blockScalarUList& result,
    const bool subtract
) const
{
    forAll(interfaces_, i)
    {
        const blockScalarList& nbr = interfaces_[i].completeExchange();
        if (subtract)
        {
            interfaces_[i].subtractCoupled(result, interfaceCoeffs_[i], nbr);
        }
        else
        {
            interfaces_[i].addCoupled(result, interfaceCoeffs_[i], nbr);
        }
    }
}


void Foam::blockLduMatrix4::Amul
(
    blockScalarUList& Ax,
    const blockScalarUList& x
) const
{
    initInterfaces(x);

    const blockScalar* __restrict__ dPtr = diag_.cdata();
    const blockScalar* __restrict__ uPtr = upper_.cdata();
    const blockScalar* __restrict__ lPtr = lower_.cdata();
    const blockScalar* __restrict__ xPtr = x.cdata();
    blockScalar* __restrict__ AxPtr = Ax.data();

    const label* const __restrict__ lAddr = addr_.lowerAddr().cdata();
    const label* const __restrict__ uAddr = addr_.upperAddr().cdata();

    for (label celli = 0; celli < nCells_; ++celli)
    {
        block4Ops::matVec
        (
            dPtr + celli*blockSize,
            xPtr + celli*blockDim,
            AxPtr + celli*blockDim
        );
    }

    for (label facei = 0; facei < nFaces_; ++facei)
    {
        const label own = lAddr[facei];
        const label nei = uAddr[facei];

        block4Ops::matVecAdd
        (
            lPtr + facei*blockSize,
            xPtr + own*blockDim,
            AxPtr + nei*blockDim
        );
        block4Ops::matVecAdd
        (
            uPtr + facei*blockSize,
            xPtr + nei*blockDim,
            AxPtr + own*blockDim
        );
    }

    updateInterfaces(Ax, false);
}


void Foam::blockLduMatrix4::initResidualAcc(const blockScalarUList& b) const
{
    const label n = nRows();
    if (residualAcc_.size() != n)
    {
        residualAcc_.resize_nocopy(n);
    }

    const blockScalar* __restrict__ bPtr = b.cdata();
    reduceScalar* __restrict__ accPtr = residualAcc_.data();
    for (label i = 0; i < n; ++i)
    {
        accPtr[i] = toDouble(bPtr[i]);
    }
}


void Foam::blockLduMatrix4::subtractAmulDouble
(
    const blockScalarUList& x
) const
{
    // Same pattern as Amul (exchange posted before the local work, completed
    // after it), but subtracting from the double accumulator
    initInterfaces(x);

    const blockScalar* __restrict__ dPtr = diag_.cdata();
    const blockScalar* __restrict__ uPtr = upper_.cdata();
    const blockScalar* __restrict__ lPtr = lower_.cdata();
    const blockScalar* __restrict__ xPtr = x.cdata();
    reduceScalar* __restrict__ accPtr = residualAcc_.data();

    const label* const __restrict__ lAddr = addr_.lowerAddr().cdata();
    const label* const __restrict__ uAddr = addr_.upperAddr().cdata();

    for (label celli = 0; celli < nCells_; ++celli)
    {
        block4Ops::matVecSubDouble
        (
            dPtr + celli*blockSize,
            xPtr + celli*blockDim,
            accPtr + celli*blockDim
        );
    }

    for (label facei = 0; facei < nFaces_; ++facei)
    {
        const label own = lAddr[facei];
        const label nei = uAddr[facei];

        block4Ops::matVecSubDouble
        (
            lPtr + facei*blockSize,
            xPtr + own*blockDim,
            accPtr + nei*blockDim
        );
        block4Ops::matVecSubDouble
        (
            uPtr + facei*blockSize,
            xPtr + nei*blockDim,
            accPtr + own*blockDim
        );
    }

    forAll(interfaces_, i)
    {
        const blockScalarList& nbr = interfaces_[i].completeExchange();
        interfaces_[i].subtractCoupledDouble
        (
            residualAcc_,
            interfaceCoeffs_[i],
            nbr
        );
    }
}


void Foam::blockLduMatrix4::narrowResidualAcc(blockScalarUList& r) const
{
    const label n = nRows();
    const reduceScalar* __restrict__ accPtr = residualAcc_.cdata();
    blockScalar* __restrict__ rPtr = r.data();
    for (label i = 0; i < n; ++i)
    {
        rPtr[i] = narrow(accPtr[i]);
    }
}


void Foam::blockLduMatrix4::residual
(
    blockScalarUList& r,
    const blockScalarUList& x,
    const blockScalarUList& b
) const
{
    Amul(r, x);

    const label n = nRows();
    blockScalar* __restrict__ rPtr = r.data();
    const blockScalar* __restrict__ bPtr = b.cdata();
    for (label i = 0; i < n; ++i)
    {
        rPtr[i] = bPtr[i] - rPtr[i];
    }
}


void Foam::blockLduMatrix4::residualDouble
(
    blockScalarUList& r,
    const UList<reduceScalar>& x,
    const blockScalarUList& b
) const
{
    const label n = nRows();
    if (residualHi_.size() != n)
    {
        residualHi_.resize_nocopy(n);
        residualLo_.resize_nocopy(n);
    }

    // x = x_hi + x_lo up to O(eps_float^2 |x|)
    {
        const reduceScalar* __restrict__ xPtr = x.cdata();
        blockScalar* __restrict__ hiPtr = residualHi_.data();
        blockScalar* __restrict__ loPtr = residualLo_.data();
        for (label i = 0; i < n; ++i)
        {
            const blockScalar hi = narrow(xPtr[i]);
            hiPtr[i] = hi;
            loPtr[i] = narrow(xPtr[i] - toDouble(hi));
        }
    }

    initResidualAcc(b);
    subtractAmulDouble(residualHi_);
    subtractAmulDouble(residualLo_);
    narrowResidualAcc(r);
}


Foam::label Foam::blockLduMatrix4::nNonFinite() const
{
    label n = 0;
    auto count = [&n](const blockScalarList& a)
    {
        for (const blockScalar v : a)
        {
            if (!std::isfinite(v)) ++n;
        }
    };
    count(diag_);
    count(upper_);
    count(lower_);
    count(source_);
    forAll(interfaceCoeffs_, i)
    {
        count(interfaceCoeffs_[i]);
    }
    return n;
}


// ************************************************************************* //
