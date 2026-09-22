/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::rhieChow

Description
    Rhie-Chow interpolation of the coupled system (spec 5.3 d, e; D-003,
    D-013).

    D_P = V_P / abar_P with abar_P the mean of the three momentum diagonal
    coefficients (incl. PTC, excl. the pressure column). Face values:
        internal       D_f = w D_P + (1-w) D_N
        coupled        same with the neighbour value across the interface
        physical       D_f = D_P where p is fixed, 0 elsewhere (D-013)

    Face flux (5.3e), identical to the continuity row of the matrix:
        phi_f = U_f.S_f - g_f (p_N - p_P) + q_f
        g_f   = D_f |S_f| deltaCoeff_f           (implicit, p-p block)
        q_f   = D_f ( gradp_f.S_f - |S_f| corr_f )   (explicit)
    with gradp_f the linear interpolate of the cell gradients ("average
    gradient") and corr_f the limited non-orthogonal correction of snGrad(p).
    On physical patches U_f, p_N are the boundary values and the explicit
    term is (1 - viC_p) D_f gradp_P.S_f with viC_p the native
    valueInternalCoeffs of p (1 for fixedValue, valueFraction for mixed
    types, 0 for zeroGradient where D_f = 0 anyway), i.e. the same weight as
    the implicit term g (1 - viC_p); corr_f = 0 there (native boundary
    snGrad). In static cells gradp is cellLimited (D-018), for the average
    gradient and for the non-orthogonal correction corr_f alike.

    Tensorial diffusivity (amendment C2, D-053), coupled.rhieChow.tensorial
    (default yes): D_P = V_P A_P^-1 with A_P the 3x3 momentum diagonal
    block (after PTC, incl. implicit MRF Coriolis, excl. the pressure
    column), in double, once per assembly; block4Ops::invert3, or the
    Moore-Penrose pseudo-inverse (block4Ops::pseudoInverse3) if
    |det A_P| < 1e-12 |trace(A_P)/3|^3, counted as nPseudoInverse. Face
    value D_f = n.interpolate(D)_f.n (n = S_f/|S_f|, linear interpolation
    of the 9 entries; neighbour values across coupled interfaces from the
    volTensorField exchange, as for the scalar D); on physical patches
    D_f = n.D_P.n where p is fixed, 0 elsewhere. D_f then enters g_f and
    q_f exactly as the scalar D_f does, so the continuity row and the flux
    update stay identical (D-035). The scalar D = V/abar is still built: it
    is the row-scaling reference Dref (5.3f) and the written coupledD.
    tensorial no: the scalar D_f of D = V/abar (bit-identical to before C2).

SourceFiles
    rhieChow.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_rhieChow_H
#define coupledFoam_rhieChow_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "boundaryCoupling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonOrthCorrection;

class rhieChow
{
    // Private Data

        const fvMesh& mesh_;

        //- Patch classification
        const List<boundaryCoupling::patchKind> kinds_;

        //- D = V/abar, cell values with coupled boundary values
        volScalarField D_;

        //- Face values D_f
        surfaceScalarField Df_;

        //- Cell gradient of p used by the last q evaluation
        volVectorField gradp_;

        //- Explicit face term q_f of the last evaluation
        surfaceScalarField q_;

        //- True once D has been set
        bool valid_;

        //- Remediation cells with gradLimiter: limited gradient of p (8.1,
        //  D-018, D-066)
        boolList isStatic_;
        bool anyStatic_;

        //- Their gradient scheme (remediation.limitedGradScheme)
        string limitedGradScheme_;

        //- Tensorial diffusivity (C2)
        bool tensorial_;
        doubleScalar detRelTol_;
        doubleScalar pinvRelTol_;
        label pinvMaxSweeps_;

        //- D = V A^-1 per cell with coupled boundary values (tensorial only)
        autoPtr<volTensorField> DTPtr_;

        //- Pseudo-inverse fallbacks (global counts): last update, since
        //  the last takePseudoInverseWindow(), whole run
        label nPinvLast_;
        label nPinvWindow_;
        label nPinvTotal_;


    // Private Member Functions

        //- D_f = n.interpolate(DT)_f.n with the boundary policy (C2)
        void buildDfTensor(const volScalarField& p);


public:

    // Constructors

        //- Construct; reads coupledDict.rhieChow (tensorial, C2)
        rhieChow(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        //- Remediation cells with gradLimiter (8.1, D-066): the pressure
        //  gradient of the Rhie-Chow terms uses remediation.limitedGradScheme
        //  (default "cellLimited Gauss linear 1") there
        void setGradLimitedCells(const boolList& isStatic);

        //- Set D from abar (per cell) and build the scalar D_f with the
        //  boundary policy (tensorial no; restart fallback)
        void updateD(const scalarField& abar, const volScalarField& p);

        //- Tensorial update (C2): scalar D from abar (row-scaling
        //  reference), D = V A^-1 from the 3x3 momentum diagonal blocks A
        //  (double), D_f = n.interpolate(D)_f.n. Collective (the
        //  pseudo-inverse count is reduced).
        void updateD
        (
            const scalarField& abar,
            const tensorField& A,
            const volScalarField& p
        );

        //- Restore the tensor D (restart consistency check) and rebuild D_f
        void setDT(const volTensorField& DT, const volScalarField& p);

        //- Tensorial diffusivity selected (coupled.rhieChow.tensorial)
        bool tensorial() const noexcept { return tensorial_; }

        //- Tensor D (tensorial only)
        const volTensorField& DT() const { return *DTPtr_; }

        //- Pseudo-inverse fallbacks of the last update (global)
        label nPseudoInverse() const noexcept { return nPinvLast_; }

        //- Pseudo-inverse fallbacks of the whole run (global)
        label nPseudoInverseTotal() const noexcept { return nPinvTotal_; }

        //- Pseudo-inverse fallbacks since the previous call (global)
        label takePseudoInverseWindow()
        {
            const label n = nPinvWindow_;
            nPinvWindow_ = 0;
            return n;
        }

        //- Set D_f directly (restart consistency check)
        bool valid() const noexcept { return valid_; }

        const volScalarField& D() const noexcept { return D_; }
        const surfaceScalarField& Df() const noexcept { return Df_; }

        //- gAverage(D) (row scaling reference, 5.3f)
        scalar Dref() const;

        //- Evaluate gradp and q_f for p (needs D_f)
        void updateExplicit
        (
            const volScalarField& p,
            const nonOrthCorrection& noc
        );

        const volVectorField& gradp() const noexcept { return gradp_; }
        const surfaceScalarField& q() const noexcept { return q_; }

        //- g_f = D_f |S_f| deltaCoeff_f (internal faces)
        scalar g(const label facei) const;

        //- g on a patch face
        scalar g(const label patchi, const label pf) const;

        //- Flux update 5.3(e) from U, p.
        //  recomputeExplicit false (outer loop): q_f of the last assembly
        //  is reused, i.e. exactly the explicit term of the continuity row
        //  that was solved, so sum_f phi_f per cell equals that row's
        //  residual at the new (U, p) - conservative to the linear-solver
        //  tolerance for omega = 1 (native: phiHbyA - pEqn.flux()).
        //  recomputeExplicit true: gradp and q are evaluated from p first
        //  (restart consistency check without a stored q).
        void updateFlux
        (
            surfaceScalarField& phi,
            const volVectorField& U,
            const volScalarField& p,
            const nonOrthCorrection& noc,
            const bool recomputeExplicit
        );

        //- Set q_f directly (restart consistency check with the q_f stored
        //  at the write time)
        void setQ(const surfaceScalarField& q) { q_ = q; }

        //- Write D_f into a field for output/restart
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
