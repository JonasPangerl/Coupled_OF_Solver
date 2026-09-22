/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockSmoother

Description
    Base class of block smoothers used on every blockGAMG level, with
    runtime selection (spec 6.1 blockGAMG.smoother:
    blockGaussSeidel | blockILU0).

SourceFiles
    blockSmoother.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockSmoother_H
#define coupledFoam_blockSmoother_H

#include "blockLduMatrix4.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockSmoother
{
protected:

    // Protected Data

        const blockLduMatrix4& matrix_;

        //- Pivot guard for the 4x4 inversions
        const reduceScalar pivotGuard_;

        //- Guarded pivots at the last update (this rank)
        label nSingular_;

        //- Under-relaxation of each sweep's update (1 = undamped)
        const blockScalar relax_;

        //- True if relax_ != 1 (the undamped path stays bit-identical)
        const bool damped_;


public:

    //- Runtime type information
    TypeName("blockSmoother");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockSmoother,
            dictionary,
            (
                const blockLduMatrix4& matrix,
                const dictionary& dict
            ),
            (matrix, dict)
        );


    // Constructors

        blockSmoother(const blockLduMatrix4& matrix, const dictionary& dict);

        //- No copy
        blockSmoother(const blockSmoother&) = delete;
        void operator=(const blockSmoother&) = delete;


    // Selectors

        //- Select from the "smoother" keyword of dict
        static autoPtr<blockSmoother> New
        (
            const blockLduMatrix4& matrix,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blockSmoother() = default;


    // Member Functions

        //- Recompute inverses/factorisation after a coefficient change
        virtual void update() = 0;

        //- nSweeps smoothing sweeps on A x = b, x updated in place
        virtual void smooth
        (
            blockScalarUList& x,
            const blockScalarUList& b,
            const label nSweeps
        ) const = 0;

        //- Guarded pivots at the last update
        label nSingularDiag() const noexcept
        {
            return nSingular_;
        }

        //- Cells whose factorised pivot fell back to the unmodified
        //  diagonal block at the last update (blockILU0 pivot-growth
        //  safeguard; 0 for the other smoothers)
        virtual label nPivotFallback() const noexcept
        {
            return 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
