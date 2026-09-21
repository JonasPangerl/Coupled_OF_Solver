/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockBiCGStab

Description
    Right-preconditioned BiCGStab (Saad 2003, Alg. 7.7; same arrangement as
    native PBiCGStab) on the 4*nCells block vector (spec 6.2).

    Breakdown guards: |rho| < VSMALL*||r0||^2 or |omega| < VSMALL restarts
    from the current x with a fresh shadow residual, at most maxRestarts
    times, then returns with converged = false and breakdown = true.

SourceFiles
    blockBiCGStab.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockBiCGStab_H
#define coupledFoam_blockBiCGStab_H

#include "blockSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockBiCGStab
:
    public blockSolver
{
public:

    //- Runtime type information
    TypeName("blockBiCGStab");


    // Constructors

        blockBiCGStab(const blockLduMatrix4& matrix, const dictionary& dict);


    //- Destructor
    virtual ~blockBiCGStab() = default;


    // Member Functions

        virtual blockSolverPerformance solve
        (
            blockScalarUList& x,
            const blockScalarUList& b,
            const reduceScalar normFactor = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
