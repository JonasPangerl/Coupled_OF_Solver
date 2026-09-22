/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGMRES

Description
    Restarted GMRES(m) with right preconditioning (Saad 2003, Alg. 9.5;
    amendment B1): Arnoldi on A M^-1, only V_j is stored and the correction
    x += M^-1 (V y) is formed once per restart cycle. Correct only for a
    fixed (linear) preconditioner; use blockFGMRES with the K-cycle.
    Happy-breakdown guard as blockFGMRES.

    Mixed-precision iterative refinement exactly as blockFGMRES (GMRES-IR):
    the Krylov basis, the preconditioner and the matrix are blockScalar, the
    iterate is kept in double and every restart cycle starts from the true
    residual b - A x of that double iterate, evaluated in double
    (blockLduMatrix4::residualDouble). Convergence is decided by that true
    residual only; finalResidual is the true residual of the double iterate
    and the returned x its rounding to blockScalar. With a blockScalar
    iterate the attainable residual was bounded by the rounding of x
    (Test-blockFGMRES, fixed preconditioner: stall at 1.6e-8).

    Dictionary: restart (default 10), legacy keyword gmresRestart.

    Memory: m + 3 block vectors of 4*nCells blockScalar, plus the double
    iterate and, inside the matrix, the double residual accumulator and two
    blockScalar split vectors.

SourceFiles
    blockGMRES.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGMRES_H
#define coupledFoam_blockGMRES_H

#include "blockSolver.H"
#include "alignedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockGMRES
:
    public blockSolver
{
    // Private Data

        //- Krylov subspace size before restart
        const label restart_;

        //- Work vectors, kept across solves: no allocation of m+4 vectors
        //  per solve (D-069 F12). Sized on first use (and when the level
        //  size changes); the contents are always written before read.
        mutable alignedList<blockScalar> r_;
        mutable alignedList<blockScalar> w_;
        mutable List<alignedList<blockScalar>> V_;
        mutable alignedList<blockScalar> z_;
        mutable alignedList<blockScalar> Vy_;
        mutable alignedList<reduceScalar> xd_;


public:

    //- Runtime type information
    TypeName("blockGMRES");


    // Constructors

        blockGMRES(const blockLduMatrix4& matrix, const dictionary& dict);


    //- Destructor
    virtual ~blockGMRES() = default;


    // Member Functions

        virtual blockSolverPerformance solve
        (
            blockScalarUList& x,
            const blockScalarUList& b,
            const reduceScalar normFactor = -1
        ) const;

        virtual void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
