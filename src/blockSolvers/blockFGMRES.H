/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockFGMRES

Description
    Flexible GMRES(m) (Saad 2003, Alg. 9.6; amendment B1): the
    preconditioned vectors z_j = M_j^-1 v_j are stored, so a preconditioner
    that changes between applications (block-GAMG K-cycle, iterative
    coarsest-level solve) is handled correctly. Required for cycleType K.

    Arnoldi with modified Gram-Schmidt; every inner product through
    doubleReduce. The Hessenberg matrix, the Givens rotations and the
    least-squares solve are in double. Happy-breakdown guard: if
    h_(j+1,j) < cfVSmall*||r0|| the cycle stops and the small system is solved.

    The first Arnoldi step of every solve also yields the preconditioner
    efficiency rho = ||v1 - A M^-1 v1|| (autoTune, 6.3.5) at the cost of one
    extra reduction.

    Mixed-precision iterative refinement (GMRES-IR, Carson & Higham 2017,
    2018): the Krylov basis, the preconditioner and the matrix are
    blockScalar, but the iterate is kept in double and every restart starts
    from the true residual b - A x of that double iterate, evaluated in
    double (blockLduMatrix4::residualDouble). Each restart cycle is a
    correction solve in blockScalar, the refinement is in double.
    Convergence is decided by that true residual only. With a blockScalar
    iterate the attainable true residual is bounded by the rounding of x
    itself, ~eps_float*|| |A| |x| || - on the Test-blockGAMG system about
    6.5e-9 in the 6.1 norm, above the 1e-9 test tolerance (D-023).
    finalResidual is the true residual of the double iterate; the returned
    x is its rounding to blockScalar.

    Dictionary: restart (default 10; 6 above 40 M cells, B7), legacy
    keyword gmresRestart accepted.

    Memory: 2*m + 2 block vectors of 4*nCells blockScalar, plus the double
    iterate and, inside the matrix, the double residual accumulator and two
    blockScalar split vectors.

SourceFiles
    blockFGMRES.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockFGMRES_H
#define coupledFoam_blockFGMRES_H

#include "blockSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockFGMRES
:
    public blockSolver
{
    // Private Data

        //- Krylov subspace size before restart
        const label restart_;


public:

    //- Runtime type information
    TypeName("blockFGMRES");


    // Constructors

        blockFGMRES(const blockLduMatrix4& matrix, const dictionary& dict);


    //- Destructor
    virtual ~blockFGMRES() = default;


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
