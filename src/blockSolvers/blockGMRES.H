/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGMRES

Description
    Restarted GMRES(m) (Saad 2003, Alg. 6.9) with right preconditioning in
    the flexible form (FGMRES, Saad Alg. 9.6): the preconditioned vectors
    z_j = M^-1 v_j are stored, so a preconditioner that is not exactly linear
    (block-GAMG with an iterative coarsest-level solve) is handled correctly.

    Arnoldi with modified Gram-Schmidt; every inner product through
    doubleReduce. The Hessenberg matrix, the Givens rotations and the
    least-squares solve are in double.

    Dictionary: gmresRestart (default 10).

    Memory: 2*m + 2 block vectors of 4*nCells blockScalar.

SourceFiles
    blockGMRES.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGMRES_H
#define coupledFoam_blockGMRES_H

#include "blockSolver.H"

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
