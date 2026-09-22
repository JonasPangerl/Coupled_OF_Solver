/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::coupledForcesFO

Description
    Function-object wrapper of coupledForces (amendment D6), so that native
    solvers (simpleFoam) write the same double-accumulated coefficients as
    coupledFoam:

        coupledForces
        {
            type            coupledForcesFO;
            libs            (coupledFoam);
            writeControl    timeStep;
            writeInterval   1;
        }

    Settings: system/coupledForcesDict (the file coupledFoam reads), or,
    if the entry itself has a `patches` keyword, the entry. Output:
    postProcessing/<entry name>/<startTime>/coeffs.dat; the iter column is
    the time value (deltaT 1: the iteration), like the Time column of the
    native forceCoeffs file.

    Inside coupledFoam the object is dormant: coupledFoam evaluates the
    same class itself after every outer iteration (no function-object
    lag) and writes the file (coupledForces::internal()).

    Results (function-object properties): Cd, Cl, CmPitch, Cd_p, Cd_v,
    Cl_p, Cl_v.

SourceFiles
    coupledForcesFO.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledForcesFO_H
#define coupledFoam_coupledForcesFO_H

#include "fvMeshFunctionObject.H"
#include "coupledForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class coupledForcesFO
:
    public fvMeshFunctionObject
{
    // Private Data

        autoPtr<coupledForces> forces_;

        //- Dormant-state message printed
        bool dormantReported_;


public:

    //- Runtime type information
    TypeName("coupledForcesFO");


    // Constructors

        coupledForcesFO
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        coupledForcesFO(const coupledForcesFO&) = delete;

        //- No copy assignment
        void operator=(const coupledForcesFO&) = delete;


    //- Destructor
    virtual ~coupledForcesFO() = default;


    // Member Functions

        virtual bool read(const dictionary& dict);

        virtual bool execute();

        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
