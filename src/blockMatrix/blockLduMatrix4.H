/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockLduMatrix4

Description
    4x4-block LDU matrix (spec 5.2) on native LDU addressing.

    Storage (all blockScalar = float):
        diag_    16*nCells   row-major 4x4 per cell, index cell*16 + r*4 + c
        upper_   16*nFaces   owner row / neighbour column (native upper)
        lower_   16*nFaces   neighbour row / owner column (native lower)
        source_   4*nCells   right-hand side
        interfaceCoeffs_[i]  16*nFaces(i) per processor interface, stored as
                             actual matrix entries (see blockLduInterface)

    Row/column meaning: 0-2 momentum x,y,z / u,v,w; 3 continuity / p.

    Amul follows native lduMatrix::Amul: diagonal, then face loop
        Ax[u[f]] += lower[f] x[l[f]]
        Ax[l[f]] += upper[f] x[u[f]]
    with the processor exchange posted before the local work and completed
    after it (native initMatrixInterfaces/updateMatrixInterfaces pattern).

    The same class is used for every GAMG level; coarse levels are built on
    the native lduPrimitiveMesh levels of GAMGAgglomeration.

SourceFiles
    blockLduMatrix4.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockLduMatrix4_H
#define coupledFoam_blockLduMatrix4_H

#include "blockScalar.H"
#include "blockLduInterface.H"
#include "lduAddressing.H"
#include "lduInterfacePtrsList.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMesh;

class blockLduMatrix4
{
    // Private Data

        //- Native addressing (owned by the mesh / agglomeration)
        const lduAddressing& addr_;

        //- The lduMesh, if constructed from one (needed by blockGAMG)
        const lduMesh* meshPtr_;

        //- Coefficient version, incremented by markUpdated()
        label version_;

        //- Communicator
        const label comm_;

        //- Sizes
        const label nCells_;
        const label nFaces_;

        //- Coefficients
        blockScalarList diag_;
        blockScalarList upper_;
        blockScalarList lower_;
        blockScalarList source_;

        //- Block-coupled (processor) interfaces
        blockLduInterfaceList interfaces_;

        //- Coefficients per block interface (same order as interfaces_)
        List<blockScalarList> interfaceCoeffs_;


public:

    // Constructors

        //- Construct zero-initialised on native addressing and interfaces.
        //  Only processor interfaces become block interfaces.
        blockLduMatrix4
        (
            const lduAddressing& addr,
            const lduInterfacePtrsList& interfaces,
            const label comm
        );

        //- Construct zero-initialised on an lduMesh
        explicit blockLduMatrix4(const lduMesh& mesh);

        //- No copy
        blockLduMatrix4(const blockLduMatrix4&) = delete;
        void operator=(const blockLduMatrix4&) = delete;


    // Member Functions

        // Access

            const lduAddressing& lduAddr() const noexcept { return addr_; }

            //- True if constructed from an lduMesh
            bool hasMesh() const noexcept { return meshPtr_ != nullptr; }

            //- The lduMesh (FatalError if constructed from addressing only)
            const lduMesh& mesh() const;

            //- Coefficient version (see markUpdated)
            label version() const noexcept { return version_; }
            label comm() const noexcept { return comm_; }
            label nCells() const noexcept { return nCells_; }
            label nFaces() const noexcept { return nFaces_; }

            //- Number of unknowns on this rank (4*nCells)
            label nRows() const noexcept { return blockDim*nCells_; }

            blockScalarList& diag() noexcept { return diag_; }
            blockScalarList& upper() noexcept { return upper_; }
            blockScalarList& lower() noexcept { return lower_; }
            blockScalarList& source() noexcept { return source_; }
            const blockScalarList& diag() const noexcept { return diag_; }
            const blockScalarList& upper() const noexcept { return upper_; }
            const blockScalarList& lower() const noexcept { return lower_; }
            const blockScalarList& source() const noexcept { return source_; }

            const blockLduInterfaceList& interfaces() const noexcept
            {
                return interfaces_;
            }

            blockScalarList& interfaceCoeffs(const label i)
            {
                return interfaceCoeffs_[i];
            }

            const blockScalarList& interfaceCoeffs(const label i) const
            {
                return interfaceCoeffs_[i];
            }

            //- Pointer to the 4x4 diagonal block of a cell
            blockScalar* diagBlock(const label celli)
            {
                return diag_.data() + celli*blockSize;
            }

            const blockScalar* diagBlock(const label celli) const
            {
                return diag_.cdata() + celli*blockSize;
            }

            //- Number of non-zero blocks on this rank
            //  (nCells + 2 nFaces + interface faces)
            label nnzBlocks() const;

            //- Bytes of coefficient storage on this rank
            std::size_t storageBytes() const;


        // Edit

            //- Zero all coefficients and the source
            void clear();

            //- Signal that the coefficients have changed (preconditioners and
            //  smoothers recompute their factorisations on the next solve)
            void markUpdated() noexcept { ++version_; }

            //- Multiply row r of every block (diag, upper, lower, interfaces)
            //  and the source row r by s. Row scaling (spec 5.3f).
            void scaleRow(const label r, const blockScalar s);


        // Operations

            //- Ax = A x (processor interfaces included)
            void Amul(blockScalarUList& Ax, const blockScalarUList& x) const;

            //- r = b - A x
            void residual
            (
                blockScalarUList& r,
                const blockScalarUList& x,
                const blockScalarUList& b
            ) const;

            //- Post the interface exchange of x (for smoothers that treat
            //  interfaces Jacobi-style)
            void initInterfaces(const blockScalarUList& x) const;

            //- Complete the exchange and add (sign +1) or subtract (sign -1)
            //  the interface contributions to result
            void updateInterfaces
            (
                blockScalarUList& result,
                const bool subtract
            ) const;

            //- Check that all stored coefficients are finite.
            //  Returns the number of non-finite entries on this rank.
            label nNonFinite() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
