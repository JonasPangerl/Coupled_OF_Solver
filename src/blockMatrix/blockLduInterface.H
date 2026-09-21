/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockLduInterface

Description
    Implicit block coupling across one processor interface (spec 5.2, 5.5
    "processor: fully implicit via block interface coefficients").

    Wraps a native lduInterface that is also a processorLduInterface
    (processorFvPatch on the finest level, processorGAMGInterface on coarse
    levels) and exchanges 4-vectors (u, v, w, p per face) as blockScalar with
    non-blocking MPI, using the native neighbour rank, tag and communicator.
    Face ordering on both sides is the native processor patch ordering.

    Sign convention: unlike native interfaceBouCoeffs (result -= coeff*psi),
    the block coefficients are stored as the actual matrix entries A_PN of
    the remote neighbour, so Amul adds coeff*x_nbr. Layout per face is the
    row-major 4x4 block, 16 values per face.

    Cyclic and cyclicAMI interfaces are not block-coupled in v2 (explicit,
    spec 5.5); isBlockCoupled() returns false for them, and also for
    processorCyclic interfaces (processor patches of a decomposed cyclic,
    lagged like cyclic in the assembly).

SourceFiles
    blockLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockLduInterface_H
#define coupledFoam_blockLduInterface_H

#include "blockScalar.H"
#include "lduInterface.H"
#include "processorLduInterface.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockLduInterface
{
    // Private Data

        //- Native interface (face cells)
        const lduInterface& interface_;

        //- Native processor interface (neighbour, tag, comm)
        const processorLduInterface& procInterface_;

        //- Index of this interface in the native interface list
        const label index_;

        //- Send buffer, 4 values per face
        mutable blockScalarList sendBuf_;

        //- Receive buffer, 4 values per face
        mutable blockScalarList recvBuf_;

        //- Outstanding requests
        mutable UPstream::Request sendReq_;
        mutable UPstream::Request recvReq_;

        //- Exchange in progress
        mutable bool pending_;


public:

    // Constructors

        //- Construct from a native interface that is a processor interface
        blockLduInterface(const lduInterface& iface, const label index);

        //- No copy
        blockLduInterface(const blockLduInterface&) = delete;
        void operator=(const blockLduInterface&) = delete;


    //- Destructor, completes any outstanding exchange
    ~blockLduInterface();


    // Static Functions

        //- True if the native interface is implicitly block-coupled
        //  (processor interfaces only, cyclic ones are explicit in v2)
        static bool isBlockCoupled(const lduInterface& iface);


    // Member Functions

        //- Native interface index
        label index() const noexcept
        {
            return index_;
        }

        //- Number of faces
        label size() const
        {
            return interface_.faceCells().size();
        }

        //- Cells adjacent to the interface faces
        const labelUList& faceCells() const
        {
            return interface_.faceCells();
        }

        //- Native interface
        const lduInterface& interface() const noexcept
        {
            return interface_;
        }

        //- Neighbour rank
        int neighbProcNo() const
        {
            return procInterface_.neighbProcNo();
        }

        //- Pack x at faceCells and post the non-blocking exchange
        void initExchange(const blockScalarUList& x) const;

        //- Wait for the exchange, return the neighbour values (4 per face)
        const blockScalarList& completeExchange() const;

        //- result[faceCells] += coeffs(16 per face) * nbr(4 per face)
        void addCoupled
        (
            blockScalarUList& result,
            const blockScalarUList& coeffs,
            const blockScalarUList& nbr
        ) const;

        //- result[faceCells] -= coeffs(16 per face) * nbr(4 per face)
        void subtractCoupled
        (
            blockScalarUList& result,
            const blockScalarUList& coeffs,
            const blockScalarUList& nbr
        ) const;

        //- As subtractCoupled, accumulated in double (residual, D-012)
        void subtractCoupledDouble
        (
            UList<reduceScalar>& result,
            const blockScalarUList& coeffs,
            const blockScalarUList& nbr
        ) const;

        //- Exchange an arbitrary per-face list of blockScalar (n per face),
        //  blocking. Used at setup time only (e.g. diagnostics).
        void exchange
        (
            const blockScalarUList& send,
            blockScalarUList& recv
        ) const;
};


//- List of block interfaces of one level
typedef PtrList<blockLduInterface> blockLduInterfaceList;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
