/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockLduInterface.H"
#include "block4Ops.H"
#include "UIPstream.H"
#include "UOPstream.H"
#include "error.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

bool Foam::blockLduInterface::isBlockCoupled(const lduInterface& iface)
{
    return dynamic_cast<const processorLduInterface*>(&iface) != nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockLduInterface::blockLduInterface
(
    const lduInterface& iface,
    const label index
)
:
    interface_(iface),
    procInterface_(dynamic_cast<const processorLduInterface&>(iface)),
    index_(index),
    sendBuf_(blockDim*iface.faceCells().size()),
    recvBuf_(blockDim*iface.faceCells().size()),
    sendReq_(),
    recvReq_(),
    pending_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockLduInterface::~blockLduInterface()
{
    if (pending_)
    {
        UPstream::waitRequest(recvReq_);
        UPstream::waitRequest(sendReq_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockLduInterface::initExchange(const blockScalarUList& x) const
{
    if (pending_)
    {
        FatalErrorInFunction
            << "Exchange already in progress on interface " << index_
            << abort(FatalError);
    }

    const labelUList& fc = faceCells();
    const label nFaces = fc.size();

    for (label i = 0; i < nFaces; ++i)
    {
        const label celli = fc[i];
        for (label k = 0; k < blockDim; ++k)
        {
            sendBuf_[i*blockDim + k] = x[celli*blockDim + k];
        }
    }

    UIPstream::read
    (
        recvReq_,
        procInterface_.neighbProcNo(),
        recvBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    UOPstream::write
    (
        sendReq_,
        procInterface_.neighbProcNo(),
        sendBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    pending_ = true;
}


const Foam::blockScalarList&
Foam::blockLduInterface::completeExchange() const
{
    if (!pending_)
    {
        FatalErrorInFunction
            << "No exchange in progress on interface " << index_
            << abort(FatalError);
    }

    UPstream::waitRequest(recvReq_);
    UPstream::waitRequest(sendReq_);
    pending_ = false;

    return recvBuf_;
}


void Foam::blockLduInterface::addCoupled
(
    blockScalarUList& result,
    const blockScalarUList& coeffs,
    const blockScalarUList& nbr
) const
{
    const labelUList& fc = faceCells();
    const label nFaces = fc.size();

    for (label i = 0; i < nFaces; ++i)
    {
        block4Ops::matVecAdd
        (
            coeffs.cdata() + i*blockSize,
            nbr.cdata() + i*blockDim,
            result.data() + fc[i]*blockDim
        );
    }
}


void Foam::blockLduInterface::subtractCoupled
(
    blockScalarUList& result,
    const blockScalarUList& coeffs,
    const blockScalarUList& nbr
) const
{
    const labelUList& fc = faceCells();
    const label nFaces = fc.size();

    for (label i = 0; i < nFaces; ++i)
    {
        block4Ops::matVecSub
        (
            coeffs.cdata() + i*blockSize,
            nbr.cdata() + i*blockDim,
            result.data() + fc[i]*blockDim
        );
    }
}


void Foam::blockLduInterface::subtractCoupledDouble
(
    UList<reduceScalar>& result,
    const blockScalarUList& coeffs,
    const blockScalarUList& nbr
) const
{
    const labelUList& fc = faceCells();
    const label nFaces = fc.size();

    for (label i = 0; i < nFaces; ++i)
    {
        block4Ops::matVecSubDouble
        (
            coeffs.cdata() + i*blockSize,
            nbr.cdata() + i*blockDim,
            result.data() + fc[i]*blockDim
        );
    }
}


void Foam::blockLduInterface::exchange
(
    const blockScalarUList& send,
    blockScalarUList& recv
) const
{
    UPstream::Request rreq;
    UPstream::Request sreq;

    UIPstream::read
    (
        rreq,
        procInterface_.neighbProcNo(),
        recv,
        procInterface_.tag(),
        procInterface_.comm()
    );

    UOPstream::write
    (
        sreq,
        procInterface_.neighbProcNo(),
        send,
        procInterface_.tag(),
        procInterface_.comm()
    );

    UPstream::waitRequest(rreq);
    UPstream::waitRequest(sreq);
}


// ************************************************************************* //
