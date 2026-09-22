/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::coupledForces

Description
    Double-accumulated force and moment coefficients (amendment D6).

    Native forceCoeffs sums p S_f and the viscous face forces in `scalar`
    (float in an SP build); with the cancellation between pressure and
    suction side that is too noisy for the 0.2 % force window of 12.3(ii).
    This class forms every face product in double from the (float or
    double) face values, sums per rank in double and reduces the sums with
    one MPI_DOUBLE call.

    Per face f of the selected patches (kinematic p, incompressible):
        F_p = rhoInf (p_f - pRef/rhoInf) S_f
        F_v = rhoInf S_f . devReff_f      (devReff = -nuEff dev(2 symm grad U)
                                           carries the minus sign, i.e.
                                           F_v = -rhoInf tau_f . S_f)
        M   = (x_f - CofR) ^ (F_p + F_v)
    which is exactly the native forces::calcForcesMoments with
    rho rhoInf (devRhoReff = -rho nuEff devTwoSymm(grad U)).

    Coefficient frame, identical to native forceCoeffs (v2606
    coordSystem::cartesian(CofR, liftDir, dragDir)):
        e3 = liftDir/|liftDir|
        e1 = dragDir with its e3 component removed, normalised
        e2 = e3 ^ e1                      (the native pitch axis)
        Cd = F.e1/(q Aref), Cl = F.e3/(q Aref), Cm = M.e2/(q Aref lRef)
        q  = 0.5 rhoInf magUInf^2
    Native v2606 ignores the legacy pitchAxis entry; it is accepted here
    for documentation and only checked against e2 (message if it differs).

    Dictionary (system/coupledForcesDict, or the function-object entry):
        patches   (body);        // wordRes, patch groups allowed
        rhoInf    1.225;
        magUInf   40;
        lRef      1.044;
        Aref      0.056;
        CofR      (0 0 0);       // mesh frame; the SP harness shifts it
                                 // together with the mesh (D5.2)
        dragDir   (1 0 0);
        liftDir   (0 0 1);
        pitchAxis (0 1 0);       // optional, see above
        pRef      0;             // optional, as native (dynamic units)
        p         p;             // optional field name

    Output (master): postProcessing/<name>/<startTime>/coeffs.dat with the
    columns  iter Cd Cl Cm Cd_p Cd_v Cl_p Cl_v  (%.10e; iter %.10g).
    An existing file of the same start time is appended to (restart from
    the same time); a restart from a later time writes a new start-time
    directory, like the native function objects.

    coupledFoam evaluates the class after every outer iteration and marks
    itself as the internal evaluator (setInternal); the function-object
    wrapper coupledForcesFO is then dormant, so that the file is written
    once. In simpleFoam the wrapper evaluates the same class.

SourceFiles
    coupledForces.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledForces_H
#define coupledFoam_coupledForces_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coupledForces
{
public:

    //- Coefficients of the last evaluation (all double, D2.2 i)
    struct coefficients
    {
        double Cd = 0;
        double Cl = 0;
        double Cm = 0;
        double Cd_p = 0;
        double Cd_v = 0;
        double Cl_p = 0;
        double Cl_v = 0;
    };


private:

    // Private Data

        const fvMesh& mesh_;

        //- Output directory name below postProcessing
        word name_;

        labelList patchIDs_;

        word pName_;

        double rhoInf_;
        double magUInf_;
        double lRef_;
        double Aref_;
        double pRef_;

        double CofR_[3];
        double e1_[3];
        double e2_[3];
        double e3_[3];

        //- Global sums of the last evaluation (mesh frame)
        double Fp_[3];
        double Fv_[3];
        double Mp_[3];
        double Mv_[3];

        coefficients c_;

        bool evaluated_;

        //- Output file (master only; empty until the first write)
        fileName file_;

        //- coupledFoam evaluates internally: the function object is dormant
        static bool internal_;


    // Private Member Functions

        void readSettings(const dictionary& dict);

        //- Open (create or append) the output file, write the header if new
        void openFile();


public:

    //- Name of the case dictionary: system/coupledForcesDict
    static const word dictName;


    // Constructors

        coupledForces
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& name = "coupledForces"
        );

        //- From system/coupledForcesDict; null if the file does not exist
        static autoPtr<coupledForces> New
        (
            const fvMesh& mesh,
            const word& name = "coupledForces"
        );

        //- True if system/coupledForcesDict exists
        static bool dictPresent(const Time& runTime);

        //- The contents of system/coupledForcesDict (fatal if missing)
        static dictionary readDict(const Time& runTime);


    // Member Functions

        //- Evaluate with the kinematic pressure and devReff (the boundary
        //  values of turbulence->devReff())
        void evaluate
        (
            const volScalarField& p,
            const volSymmTensorField& devReff
        );

        //- Evaluate with the pressure field named in the dictionary
        void evaluate(const volSymmTensorField& devReff);

        //- Append one row (master): iter, Cd, Cl, Cm, Cd_p, Cd_v, Cl_p, Cl_v
        void write(const double iter);

        const coefficients& coeffs() const noexcept { return c_; }
        double Cd() const noexcept { return c_.Cd; }
        double Cl() const noexcept { return c_.Cl; }
        double Cm() const noexcept { return c_.Cm; }
        bool evaluated() const noexcept { return evaluated_; }

        const word& pName() const noexcept { return pName_; }
        const word& name() const noexcept { return name_; }

        //- Force (pressure, viscous) and moment sums, global, mesh frame
        vector pressureForce() const;
        vector viscousForce() const;
        vector pressureMoment() const;
        vector viscousMoment() const;

        //- Settings (for the effective-settings print and summaries)
        void writeSettings(dictionary& dict) const;

        //- Internal-evaluator flag (coupledFoam sets it)
        static void setInternal(const bool on) noexcept { internal_ = on; }
        static bool internal() noexcept { return internal_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
