/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGAMGPrecon

Description
    Block-GAMG preconditioner: one V-cycle of blockGAMG per application
    (spec 6.3 "V-cycle as preconditioner, not as a standalone solver").
    Selected with "preconditioner blockGAMG;"; settings from the solver
    dict's "blockGAMG" sub-dictionary (spec 6.1).

SourceFiles
    blockGAMGPrecon.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGAMGPrecon_H
#define coupledFoam_blockGAMGPrecon_H

#include "blockPreconditioner.H"
#include "blockGAMG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockGAMGPrecon
:
    public blockPreconditioner
{
    // Private Data

        //- The hierarchy
        blockGAMG gamg_;


public:

    //- Runtime type information (selection name "blockGAMG")
    TypeName("blockGAMG");


    // Constructors

        blockGAMGPrecon(const blockSolver& solver, const dictionary& dict);


    //- Destructor
    virtual ~blockGAMGPrecon() = default;


    // Member Functions

        //- The hierarchy
        const blockGAMG& gamg() const noexcept
        {
            return gamg_;
        }

        virtual void update();

        virtual void precondition
        (
            blockScalarUList& w,
            const blockScalarUList& r
        ) const;

        virtual label nSingularDiag() const
        {
            return gamg_.nSingularDiag();
        }

        virtual void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
