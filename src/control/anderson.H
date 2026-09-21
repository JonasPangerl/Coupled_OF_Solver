/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::anderson

Description
    Anderson acceleration of the outer iteration (spec amendment B5,
    Section 7.5), Walker & Ni (2011) Type II form, default off.

    The outer iteration is read as the fixed-point map
        x -> G(x) = x + f(x),   f(x) = omega*dx(x)
    where x = (U, p) is the state at which the system was assembled and
    f = omega*dx is the accepted increment after line search and clipping.
    With the stored differences of consecutive iterates and increments
        dX_j = x_(j+1) - x_j,   dF_j = f_(j+1) - f_j   (j = k-n .. k-1)
    the coefficients solve the least-squares problem
        gamma = argmin || S (f_k - sum_j gamma_j dF_j) ||_2
    and the accelerated state is
        x_(k+1) = x_k + beta*f_k - sum_j gamma_j (dX_j + beta*dF_j).
    For beta = 1 this is x_(k+1) = G(x_k) - sum_j gamma_j dG_j.

    Note on the spec text: B5 writes x_(k+1) = x_k + beta*(w_k - sum a_j dw_j)
    with x_k the state AFTER the update. That short form omits the state
    differences dX_j and, read literally with the post-update state, would
    count w_k twice. The Type II update above is the mathematically correct
    one and is what is implemented. In terms of the post-update state
    y_k = x_k + f_k that the caller passes in:
        x_(k+1) = y_k - (1 - beta)*f_k - sum_j gamma_j (dX_j + beta*dF_j).

    Norm: S scales the velocity components by 1/Uref and the pressure by
    1/pref, so that the least-squares problem is dimensionless. All dot
    products are accumulated in double and summed over all ranks.

    Least squares: thin QR of the scaled difference matrix S*dF, kept up to
    date column by column (Walker & Ni 2011, Section 5): a new column is
    orthogonalised by classical Gram-Schmidt with one re-orthogonalisation
    pass (CGS2, two global reductions); the oldest column is removed with
    Givens rotations. gamma follows from R gamma = Q^T S f_k by back
    substitution. Only Q (not S*dF) and D_j = dX_j + beta*dF_j are stored.
    A new column whose orthogonal part is below sqrt(eps_double) of its
    norm is numerically dependent (numerical-rank criterion, not a tuning
    parameter) and is not added.

    Safeguards (B5, all mandatory):
      - status skipped: max_j |gamma_j| > maxAlpha, the update is not
        applied, the history is kept (not flushed);
      - the caller checks the extrapolated state with the sentinel (9.3);
        on failure it calls reject(), which restores the un-extrapolated
        state y_k and flushes the history;
      - the caller calls flush() on every CFL cut, every change of the
        remediation set and every rollback.

    Memory (enabled only): (2*m + 2) vectors of 4 doubles per cell
    (m columns of Q, m columns of D, previous iterate and increment),
    i.e. 320 B/cell for m = 4. Nothing is allocated when disabled.

    Dictionary: coupled.anderson { enabled m beta maxAlpha }.

    References:
      Walker, H.F., Ni, P., Anderson acceleration for fixed-point
      iterations, SIAM J. Numer. Anal. 49 (2011) 1715-1735.
      Pollock, S., Rebholz, L.G., Xiao, M., Anderson-accelerated
      convergence of Picard iterations for incompressible Navier-Stokes
      equations, SIAM J. Numer. Anal. 57 (2019) 615-637.

Integration
    In applications/coupledFoam/coupledFoam.C (outer loop), with
    "anderson aa(mesh, coupledDict);" constructed next to "sentinel sen":

    1. After the assemble/solve/line-search loop: if cuts > 0 or skipStep
       (the CFL was cut), call aa.flush(). Also call aa.flush() when the
       remediation set changed in this iteration (e.g. rem.markDynamic() in
       the omegaMin branch; compare rem.nStatic()+rem.nDynamic() before and
       after rem.updateDynamic() of the previous iteration).

    2. In the field update block, while adding omega*dx to U and p, fill
         vectorField dU(mesh.nCells());  scalarField dp(mesh.nCells());
         dU[celli] = omega*vector(d[0], d[1], d[2]);
         dp[celli] = omega*scalar(d[blockP]);
       (the increment actually applied, after rem.clipIncrement).

    3. After U.correctBoundaryConditions()/p.correctBoundaryConditions()
       and BEFORE assembler.rc().updateFlux(...) (so that the flux is
       computed once from the final state), if !skipStep:
         const anderson::status st =
             aa.apply(U, p, dU, dp, ls.Uref(), ls.pref());
         if (st == anderson::status::applied)
         {
             const sentinel::checkResult chkA =
                 sen.check(U, p, kPtr, omegaPtr, ls.Uref(), ls.pref());
             if (!chkA.ok)
             {
                 aa.reject(U, p);    // restores y_k, flushes
             }
         }
       If the call is placed after updateFlux instead, the flux must be
       recomputed after apply() and after reject().
       sen.store() keeps the pre-update state, so a later sentinel rollback
       still returns to x_k.

    4. In the rollback lambda (sentinel failure after the update or after
       the turbulence correction): call aa.flush().

    5. Log: aa.lastStatus(), aa.lastMaxAlpha(), aa.nApplied(),
       aa.nSkipped(), aa.nRejected(); settings: aa.writeSettings(eff).

SourceFiles
    anderson.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_anderson_H
#define coupledFoam_anderson_H

#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class anderson
{
public:

    //- Result of apply()
    enum class status
    {
        applied = 0,        //!< U, p replaced by the accelerated state
        skipped = 1,        //!< max|gamma| > maxAlpha, history kept
        collecting = 2,     //!< not enough history yet, pair stored
        disabled = 3        //!< anderson.enabled no
    };


private:

    // Private Data

        const fvMesh& mesh_;

        bool enabled_;
        label m_;
        doubleScalar beta_;
        doubleScalar maxAlpha_;

        //- Number of doubles per vector (4 per cell: u, v, w, p)
        label n_;

        // History

            //- Orthonormal columns of the thin QR of S*dF
            List<List<doubleScalar>> Q_;

            //- Update columns D_j = dX_j + beta*dF_j (unscaled)
            List<List<doubleScalar>> D_;

            //- Upper triangular R, row-major m x m
            List<doubleScalar> R_;

            //- Number of valid columns
            label nHist_;

            //- Previous iterate x_k and increment f_k
            List<doubleScalar> xPrev_;
            List<doubleScalar> fPrev_;
            bool hasPrev_;

            //- Scaling of the history (1/Uref, 1/pref)
            doubleScalar sU_;
            doubleScalar sP_;

        // Statistics

            label nApplied_;
            label nSkipped_;
            label nRejected_;
            label nFlushed_;
            status lastStatus_;
            doubleScalar lastMaxAlpha_;


    // Private Member Functions

        doubleScalar& R(const label i, const label j)
        {
            return R_[i*m_ + j];
        }

        doubleScalar R(const label i, const label j) const
        {
            return R_[i*m_ + j];
        }

        //- Allocate the history storage (first use)
        void allocate();

        //- Remove the oldest column (Givens downdate of the QR)
        void removeOldest();

        //- Orthogonalise column c of Q against columns 0..c-1 (CGS2),
        //  fill R(0..c, c). Returns false if numerically dependent.
        bool orthogonalise(const label c);

        //- Sum n doubles over all ranks in one reduction
        static void sumAll(List<doubleScalar>& v);


public:

    // Constructors

        anderson(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        bool enabled() const noexcept { return enabled_; }
        label m() const noexcept { return m_; }
        label nHistory() const noexcept { return nHist_; }
        label nApplied() const noexcept { return nApplied_; }
        label nSkipped() const noexcept { return nSkipped_; }
        label nRejected() const noexcept { return nRejected_; }
        label nFlushed() const noexcept { return nFlushed_; }
        status lastStatus() const noexcept { return lastStatus_; }
        doubleScalar lastMaxAlpha() const noexcept { return lastMaxAlpha_; }

        //- Name of a status for the log
        static word statusName(const status s);

        //- Accelerate the accepted outer update.
        //  U, p: state y_k after the update (U = x_k + dU, p = x_k + dp);
        //  dU, dp: the accepted increment f_k = omega*dx (per cell);
        //  Uref, pref: scales of the least-squares norm.
        //  On status applied, U and p hold x_(k+1) (internal field) and
        //  their boundary conditions are corrected; the flux is not.
        status apply
        (
            volVectorField& U,
            volScalarField& p,
            const vectorField& dU,
            const scalarField& dp,
            const scalar Uref,
            const scalar pref
        );

        //- The extrapolated state failed the sentinel: restore the
        //  un-extrapolated state y_k into U, p and flush the history
        void reject(volVectorField& U, volScalarField& p);

        //- Clear the history (CFL cut, remediation change, rollback)
        void flush();

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
