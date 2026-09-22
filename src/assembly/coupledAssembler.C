/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledAssembler.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "doubleReduce.H"
#include "fvCFD.H"
#include "laplacianScheme.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "findRefCell.H"
#include "IStringStream.H"
#include <limits>

using namespace Foam::boundaryCoupling;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::coupledAssembler::zeroFluxPatch(const fvPatch& patch)
{
    return
        isA<symmetryFvPatch>(patch)
     || isA<symmetryPlaneFvPatch>(patch);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledAssembler::coupledAssembler
(
    const fvMesh& mesh,
    const dictionary& coupledDict,
    const volScalarField& p,
    const MRFCoupling* mrfPtr
)
:
    timing_(false),
    times_(),
    mesh_(mesh),
    kinds_(boundaryCoupling::kinds(mesh)),
    blockIndex_(mesh.boundary().size(), -1),
    conv_(mesh),
    noc_
    (
        mesh,
        coupledDict.getOrDefault<scalar>
        (
            "nonOrthLimiter",
            coupledDefaults::nonOrthLimiter
        ),
        coupledDict.subOrEmptyDict("remediation").subOrEmptyDict("static")
            .getOrDefault<scalar>
            (
                "nonOrthLimiter",
                coupledDefaults::staticNonOrthLimiter
            )
    ),
    rc_(mesh, coupledDict),
    mrfPtr_(mrfPtr),
    sfdChi_(0),
    sfdUbarPtr_(nullptr),
    clampValue_
    (
        coupledDict.subOrEmptyDict("guards").getOrDefault<doubleScalar>
        (
            "clampValue",
            coupledDefaults::clampValue
        )
    ),
    needRef_(false),
    pRefCell_(-1),
    pRefValue_(0),
    refFluxChecked_(false),
    A_(mesh),
    rhs_(blockDim*mesh.nCells(), Zero),
    DdStage_(zeroCopy ? 0 : blockSize*mesh.nCells(), Zero),
    Ax_(blockDim*mesh.nCells(), Zero),
    b_(blockDim*mesh.nCells(), Zero),
    rMom_(mesh.nCells(), Zero),
    aMom_(mesh.nCells(), Zero),
    abar_(mesh.nCells(), Zero),
    UPtr_(nullptr),
    pPtr_(nullptr),
    sp_(1),
    normFactor_(1),
    residualL2_(0),
    rU_(0),
    rp_(0),
    nClamped_(0)
{
    // GUARD: the clamp must be representable in blockScalar, otherwise a
    // clamped value would still overflow on narrowing (9.2)
    if
    (
        !(clampValue_ > 0)
     || clampValue_ > doubleScalar(std::numeric_limits<blockScalar>::max())
    )
    {
        FatalIOErrorInFunction(coupledDict)
            << "guards.clampValue " << clampValue_ << " must be in (0, "
            << std::numeric_limits<blockScalar>::max()
            << "] (largest blockScalar)" << exit(FatalIOError);
    }

    forAll(A_.interfaces(), bi)
    {
        blockIndex_[A_.interfaces()[bi].index()] = bi;
    }

    forAll(kinds_, patchi)
    {
        if (kinds_[patchi] == patchKind::processor && blockIndex_[patchi] < 0)
        {
            FatalErrorInFunction
                << "Processor patch " << mesh.boundary()[patchi].name()
                << " has no block interface" << exit(FatalError);
        }
    }

    // Pressure reference for closed domains (native setRefCell)
    if (p.needReference())
    {
        needRef_ = true;
        setRefCell(p, coupledDict, pRefCell_, pRefValue_);
    }
}


void Foam::coupledAssembler::checkReferenceFluxBalance
(
    const volVectorField& U
)
{
    // adjustPhi-style check (native adjustPhi): with the strong pressure
    // reference (D-021) one continuity equation is dropped, which is exact
    // only if the prescribed boundary fluxes sum to zero. Collective:
    // needRef_ is the same on all ranks (p.needReference() reduces).
    reduceScalar sums[2] = {0, 0};   // sum U_b.S_f, sum |U_b.S_f|
    const surfaceVectorField& Sf = mesh_.Sf();
    forAll(kinds_, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        if (kinds_[patchi] != patchKind::physical || zeroFluxPatch(patch))
        {
            continue;
        }
        const vectorField& Ub = U.boundaryField()[patchi];
        const vectorField& Sp = Sf.boundaryField()[patchi];
        forAll(Ub, pf)
        {
            const reduceScalar f = Ub[pf] & Sp[pf];
            sums[0] += f;
            sums[1] += std::abs(f);
        }
    }
    doubleReduce::parSum(sums, 2, mesh_.comm());

    if
    (
        std::abs(sums[0])
      > coupledDefaults::refFluxBalanceTol*max(sums[1], cfVSmall<reduceScalar>())
    )
    {
        WarningInFunction
            << "Closed domain (pressure reference, D-021) but the boundary"
            << " fluxes do not balance: sum U_b.S_f = " << sums[0]
            << ", sum |U_b.S_f| = " << sums[1] << " (relative tolerance "
            << coupledDefaults::refFluxBalanceTol << "). The continuity"
            << " equation of the reference cell is dropped, so the"
            << " imbalance is silently absorbed there." << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledAssembler::setStaticCells(const boolList& isStatic)
{
    noc_.setStaticCells(isStatic);
    rc_.setStaticCells(isStatic);
}


void Foam::coupledAssembler::assembleMomentum
(
    const volVectorField& U,
    const volScalarField& p,
    const surfaceScalarField& phi,
    const volScalarField& nuEff,
    const scalarField& beta
)
{
    const doubleScalar tStart = (timing_ ? diagnostics::clock() : 0);

    UPtr_ = &U;
    pPtr_ = &p;

    A_.clear();
    DdStage_ = Zero;
    diagScalar* __restrict__ Dd = diagAcc();
    Ax_ = Zero;
    b_ = Zero;
    nClamped_ = 0;

    const scalarField& V = mesh_.V();
    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();

    // --- Implicit momentum operator: upwind convection - uncorrected
    //     Laplacian (non-orthogonal part explicit below, D-014)
    const surfaceScalarField nuEfff("interpolate(nuEff)", fvc::interpolate(nuEff));

    tmp<fvVectorMatrix> tUEqn = conv_.implicitPart(phi, U);
    {
        IStringStream schemeData("Gauss linear uncorrected");
        tmp<fv::laplacianScheme<vector, scalar>> tlap =
            fv::laplacianScheme<vector, scalar>::New(mesh_, schemeData);
        tUEqn.ref() -= tlap.ref().fvmLaplacian(nuEfff, U);
    }
    const fvVectorMatrix& UEqn = tUEqn();

    // --- Explicit momentum sources per unit volume (right-hand side sign)
    vectorField Sexp(-conv_.deferredCorrection(phi, U, beta));
    {
        // Transpose stress, as native linearViscousStress
        const volTensorField tau
        (
            "(nuEff*dev2(T(grad(U))))",
            nuEff*dev2(T(fvc::grad(U)))
        );
        Sexp += fvc::div(tau)().primitiveField();
    }
    if (!noc_.orthogonal())
    {
        // Limited non-orthogonal correction of the Laplacian; physical
        // boundary faces carry none (native boundary snGrad)
        surfaceVectorField corrFlux
        (
            "coupledAssembler::corrFlux",
            nuEfff*mesh_.magSf()*noc_.correction(U)
        );
        forAll(kinds_, patchi)
        {
            if (kinds_[patchi] == patchKind::physical)
            {
                corrFlux.boundaryFieldRef()[patchi] = Zero;
            }
        }
        Sexp += fvc::div(corrFlux)().primitiveField();
    }

    // --- Cells: momentum diagonal, source, MRF
    const scalarField& diagU = UEqn.diag();
    const vectorField& srcU = UEqn.source();

    forAll(V, celli)
    {
        for (label c = 0; c < blockP; ++c)
        {
            Dd[di(celli, c, c)] += diagU[celli];
            b_[celli*blockDim + c] += srcU[celli][c] + Sexp[celli][c]*V[celli];
        }
    }

    if (mrfPtr_ && mrfPtr_->active())
    {
        const tensorField& C = mrfPtr_->coeffs();
        forAll(C, celli)
        {
            for (label r = 0; r < blockP; ++r)
            {
                for (label c = 0; c < blockP; ++c)
                {
                    Dd[di(celli, r, c)] += C[celli][r*blockP + c];
                }
            }
        }
    }

    // --- SFD forcing -chi*(U - Ubar) (7.6), before PTC
    if (sfdChi_ > 0 && sfdUbarPtr_)
    {
        const vectorField& Ub = *sfdUbarPtr_;
        forAll(V, celli)
        {
            const doubleScalar cv = doubleScalar(sfdChi_)*V[celli];
            for (label c = 0; c < blockP; ++c)
            {
                Dd[di(celli, c, c)] += diagScalar(cv);
                b_[celli*blockDim + c] += cv*Ub[celli][c];
            }
        }
    }

    // --- Internal faces: momentum coupling and pressure gradient
    const scalarField& upperU = UEqn.upper();
    const scalarField& lowerU = UEqn.lower();

    blockScalar* __restrict__ Aup = A_.upper().data();
    blockScalar* __restrict__ Alo = A_.lower().data();

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wf = w[facei];
        const vector& S = Sf[facei];
        blockScalar* Ub = Aup + facei*blockSize;
        blockScalar* Lb = Alo + facei*blockSize;

        for (label c = 0; c < blockP; ++c)
        {
            const doubleScalar uc = upperU[facei];
            const doubleScalar lc = lowerU[facei];
            Ub[c*blockDim + c] = store(uc);
            Lb[c*blockDim + c] = store(lc);
            Ax_[P*blockDim + c] += uc*Ui[N][c];
            Ax_[N*blockDim + c] += lc*Ui[P][c];

            // Green-Gauss pressure gradient, p_f = w p_P + (1-w) p_N
            Dd[di(P, c, blockP)] += wf*S[c];
            const doubleScalar ug = (1 - wf)*S[c];
            Ub[c*blockDim + blockP] = store(ug);
            Ax_[P*blockDim + c] += ug*pi[N];

            Dd[di(N, c, blockP)] -= (1 - wf)*S[c];
            const doubleScalar lg = -wf*S[c];
            Lb[c*blockDim + blockP] = store(lg);
            Ax_[N*blockDim + c] += lg*pi[P];
        }
    }

    // --- Boundary faces
    const doubleScalar tb0 = (timing_ ? diagnostics::clock() : 0);
    forAll(kinds_, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& fc = patch.faceCells();
        const vectorField& Sp = Sf.boundaryField()[patchi];
        const scalarField& wp = w.boundaryField()[patchi];
        const vectorField& ic = UEqn.internalCoeffs()[patchi];
        const vectorField& bc = UEqn.boundaryCoeffs()[patchi];

        switch (kinds_[patchi])
        {
            case patchKind::empty:
                break;

            case patchKind::processor:
            {
                blockScalarList& C = A_.interfaceCoeffs(blockIndex_[patchi]);
                const vectorField Un
                (
                    U.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField pn
                (
                    p.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(fc, pf)
                {
                    const label P = fc[pf];
                    blockScalar* Cb = C.data() + pf*blockSize;
                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, c, c)] += ic[pf][c];
                        // Native sign: result -= bouCoeffs*psi_nbr
                        const doubleScalar cc = -bc[pf][c];
                        Cb[c*blockDim + c] = store(cc);
                        Ax_[P*blockDim + c] += cc*Un[pf][c];

                        Dd[di(P, c, blockP)] += wp[pf]*Sp[pf][c];
                        const doubleScalar gc = (1 - wp[pf])*Sp[pf][c];
                        Cb[c*blockDim + blockP] = store(gc);
                        Ax_[P*blockDim + c] += gc*pn[pf];
                    }
                }
                break;
            }

            case patchKind::explicitCoupled:
            {
                const vectorField Un
                (
                    U.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField pn
                (
                    p.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(fc, pf)
                {
                    const label P = fc[pf];
                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, c, c)] += ic[pf][c];
                        b_[P*blockDim + c] += bc[pf][c]*Un[pf][c];

                        Dd[di(P, c, blockP)] += wp[pf]*Sp[pf][c];
                        b_[P*blockDim + c] -= (1 - wp[pf])*Sp[pf][c]*pn[pf];
                    }
                }
                break;
            }

            case patchKind::physical:
            {
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const scalarField viCp(pp.valueInternalCoeffs(wp));
                const scalarField vbCp(pp.valueBoundaryCoeffs(wp));

                forAll(fc, pf)
                {
                    const label P = fc[pf];
                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, c, c)] += ic[pf][c];
                        b_[P*blockDim + c] += bc[pf][c];

                        // p_f = viC p_P + vbC
                        Dd[di(P, c, blockP)] += Sp[pf][c]*viCp[pf];
                        b_[P*blockDim + c] -= Sp[pf][c]*vbCp[pf];
                    }
                }
                break;
            }
        }
    }

    const doubleScalar tb1 = (timing_ ? diagnostics::clock() : 0);

    // --- Diagonal contribution to A x for rows 0-2 (PTC excluded: it
    //     cancels in b - A x), momentum residual
    forAll(V, celli)
    {
        const doubleScalar x[blockDim] =
            {Ui[celli][0], Ui[celli][1], Ui[celli][2], pi[celli]};

        doubleScalar d = 0;
        for (label r = 0; r < blockP; ++r)
        {
            doubleScalar s = 0;
            for (label c = 0; c < blockDim; ++c)
            {
                s += Dd[di(celli, r, c)]*x[c];
            }
            Ax_[celli*blockDim + r] += s;
            rMom_[celli][r] =
                scalar(b_[celli*blockDim + r] - Ax_[celli*blockDim + r]);
            d += Dd[di(celli, r, r)];
        }
        aMom_[celli] = scalar(d/blockP);
    }

    if (timing_)
    {
        const doubleScalar tEnd = diagnostics::clock();
        times_.boundary += tb1 - tb0;
        times_.momentumOps += (tEnd - tStart) - (tb1 - tb0);
    }
}


void Foam::coupledAssembler::assembleContinuity(const scalarField& rDeltaTV)
{
    const doubleScalar tStart = (timing_ ? diagnostics::clock() : 0);

    if (!UPtr_ || !pPtr_)
    {
        FatalErrorInFunction
            << "assembleMomentum must be called first" << abort(FatalError);
    }

    const volVectorField& U = *UPtr_;
    const volScalarField& p = *pPtr_;

    const surfaceScalarField& w = mesh_.weights();
    const surfaceVectorField& Sf = mesh_.Sf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Ui = U.primitiveField();
    const scalarField& pi = p.primitiveField();
    const label nCells = mesh_.nCells();
    diagScalar* __restrict__ Dd = diagAcc();

    if (needRef_ && !refFluxChecked_)
    {
        checkReferenceFluxBalance(U);
        refFluxChecked_ = true;
    }

    // --- PTC diagonal and abar
    for (label celli = 0; celli < nCells; ++celli)
    {
        doubleScalar s = 0;
        for (label c = 0; c < blockP; ++c)
        {
            Dd[di(celli, c, c)] += rDeltaTV[celli];
            s += Dd[di(celli, c, c)];
        }
        // The SFD term is not part of abar: the Rhie-Chow D = V/abar of a
        // converged state must not depend on chi (D-052)
        abar_[celli] = scalar
        (
            s/blockP
          - (sfdUbarPtr_ ? doubleScalar(sfdChi_)*mesh_.V()[celli] : 0)
        );
    }

    // --- Rhie-Chow D, D_f, explicit face term q
    const doubleScalar tr0 = (timing_ ? diagnostics::clock() : 0);
    if (rc_.tensorial())
    {
        // C2: 3x3 momentum diagonal block after PTC (incl. MRF Coriolis,
        // excl. the pressure column), double
        // (C3: without the SFD term, D-052)
        tensorField Amom(nCells);
        for (label celli = 0; celli < nCells; ++celli)
        {
            tensor& a = Amom[celli];
            for (label r = 0; r < blockP; ++r)
            {
                for (label c = 0; c < blockP; ++c)
                {
                    a[r*blockP + c] = Dd[di(celli, r, c)];
                }
            }
            if (sfdUbarPtr_)
            {
                const doubleScalar cv =
                    doubleScalar(sfdChi_)*mesh_.V()[celli];
                for (label r = 0; r < blockP; ++r)
                {
                    a[r*blockP + r] -= scalar(cv);
                }
            }
        }
        rc_.updateD(abar_, Amom, p);
    }
    else
    {
        rc_.updateD(abar_, p);
    }
    rc_.updateExplicit(p, noc_);
    const doubleScalar tr1 = (timing_ ? diagnostics::clock() : 0);

    // --- Continuity row scale s_p (5.3f); GUARD: s_p <= 1/cfVSmall. Every
    //     row-3 coefficient is scaled in double before it is narrowed, so
    //     the clamp and the single rounding apply to the scaled value.
    sp_ = 1.0/max(rc_.Dref(), cfVSmall<doubleScalar>());
    const doubleScalar sp = sp_;
    const surfaceScalarField& q = rc_.q();

    blockScalar* __restrict__ Aup = A_.upper().data();
    blockScalar* __restrict__ Alo = A_.lower().data();

    // --- Internal faces, row 3: sum_f phi_f = 0 with
    //     phi_f = w S.U_P + (1-w) S.U_N - g (p_N - p_P) + q
    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wf = w[facei];
        const vector& S = Sf[facei];
        const doubleScalar g = rc_.g(facei);
        blockScalar* Ub = Aup + facei*blockSize;
        blockScalar* Lb = Alo + facei*blockSize;

        for (label c = 0; c < blockP; ++c)
        {
            // Owner row (+)
            Dd[di(P, blockP, c)] += wf*S[c];
            const doubleScalar ud = (1 - wf)*S[c];
            Ub[blockP*blockDim + c] = store(sp*ud);
            Ax_[P*blockDim + blockP] += ud*Ui[N][c];

            // Neighbour row (-)
            Dd[di(N, blockP, c)] -= (1 - wf)*S[c];
            const doubleScalar ld = -wf*S[c];
            Lb[blockP*blockDim + c] = store(sp*ld);
            Ax_[N*blockDim + blockP] += ld*Ui[P][c];
        }

        Dd[di(P, blockP, blockP)] += diagScalar(g);
        Ub[blockP*blockDim + blockP] = store(-sp*g);
        Ax_[P*blockDim + blockP] -= g*pi[N];

        Dd[di(N, blockP, blockP)] += diagScalar(g);
        Lb[blockP*blockDim + blockP] = store(-sp*g);
        Ax_[N*blockDim + blockP] -= g*pi[P];

        b_[P*blockDim + blockP] -= q[facei];
        b_[N*blockDim + blockP] += q[facei];
    }

    // --- Boundary faces, row 3
    const doubleScalar tb0 = (timing_ ? diagnostics::clock() : 0);
    forAll(kinds_, patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& fc = patch.faceCells();
        const vectorField& Sp = Sf.boundaryField()[patchi];
        const scalarField& wp = w.boundaryField()[patchi];
        const scalarField& qp = q.boundaryField()[patchi];

        switch (kinds_[patchi])
        {
            case patchKind::empty:
                break;

            case patchKind::processor:
            {
                blockScalarList& C = A_.interfaceCoeffs(blockIndex_[patchi]);
                const vectorField Un
                (
                    U.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField pn
                (
                    p.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(fc, pf)
                {
                    const label P = fc[pf];
                    const doubleScalar g = rc_.g(patchi, pf);
                    blockScalar* Cb = C.data() + pf*blockSize;

                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, blockP, c)] += wp[pf]*Sp[pf][c];
                        const doubleScalar cd = (1 - wp[pf])*Sp[pf][c];
                        Cb[blockP*blockDim + c] = store(sp*cd);
                        Ax_[P*blockDim + blockP] += cd*Un[pf][c];
                    }

                    Dd[di(P, blockP, blockP)] += diagScalar(g);
                    Cb[blockP*blockDim + blockP] = store(-sp*g);
                    Ax_[P*blockDim + blockP] -= g*pn[pf];

                    b_[P*blockDim + blockP] -= qp[pf];
                }
                break;
            }

            case patchKind::explicitCoupled:
            {
                const vectorField Un
                (
                    U.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField pn
                (
                    p.boundaryField()[patchi].patchNeighbourField()
                );

                forAll(fc, pf)
                {
                    const label P = fc[pf];
                    const doubleScalar g = rc_.g(patchi, pf);

                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, blockP, c)] += wp[pf]*Sp[pf][c];
                    }
                    b_[P*blockDim + blockP] -=
                        (1 - wp[pf])*(Sp[pf] & Un[pf]);

                    Dd[di(P, blockP, blockP)] += diagScalar(g);
                    b_[P*blockDim + blockP] += g*pn[pf];

                    b_[P*blockDim + blockP] -= qp[pf];
                }
                break;
            }

            case patchKind::physical:
            {
                if (zeroFluxPatch(patch))
                {
                    // U_f.S_f = 0 exactly, no Rhie-Chow term (D_f = 0)
                    break;
                }

                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const vectorField viCU(Up.valueInternalCoeffs(wp));
                const vectorField vbCU(Up.valueBoundaryCoeffs(wp));
                const scalarField viCp(pp.valueInternalCoeffs(wp));
                const scalarField vbCp(pp.valueBoundaryCoeffs(wp));

                forAll(fc, pf)
                {
                    const label P = fc[pf];

                    // U_f = viC . U_P + vbC
                    for (label c = 0; c < blockP; ++c)
                    {
                        Dd[di(P, blockP, c)] += Sp[pf][c]*viCU[pf][c];
                    }
                    b_[P*blockDim + blockP] -= Sp[pf] & vbCU[pf];

                    // -g (p_b - p_P), p_b = viC p_P + vbC (g = 0 unless
                    // p is fixed on this patch, D-013)
                    const doubleScalar g = rc_.g(patchi, pf);
                    Dd[di(P, blockP, blockP)] += diagScalar(g*(1 - viCp[pf]));
                    b_[P*blockDim + blockP] += g*vbCp[pf];

                    b_[P*blockDim + blockP] -= qp[pf];
                }
                break;
            }
        }
    }

    const doubleScalar tb1 = (timing_ ? diagnostics::clock() : 0);

    // --- Pressure reference for closed domains (DECISIONS.md D-021):
    //     the continuity row of the reference cell is replaced by
    //     d (p_ref - pRefValue) = 0. Exact, because in a closed domain the
    //     continuity equations sum to zero and one of them is redundant;
    //     it removes the near-null space that the native weak reference
    //     (diagonal doubling) leaves for a float Krylov solve.
    if (needRef_ && pRefCell_ >= 0)
    {
        const label c = pRefCell_;
        const doubleScalar d = Dd[di(c, blockP, blockP)];

        for (label k = 0; k < blockDim; ++k)
        {
            Dd[di(c, blockP, k)] = 0;
        }
        Dd[di(c, blockP, blockP)] = diagScalar(d);

        const label nInternal = mesh_.nInternalFaces();
        for (const label facei : mesh_.cells()[c])
        {
            if (facei < nInternal)
            {
                blockScalar* blk =
                    (own[facei] == c ? Aup : Alo) + facei*blockSize;
                for (label k = 0; k < blockDim; ++k)
                {
                    blk[blockP*blockDim + k] = 0;
                }
            }
            else
            {
                const label patchi =
                    mesh_.boundaryMesh().whichPatch(facei);
                if (kinds_[patchi] == patchKind::processor)
                {
                    const label pf =
                        facei - mesh_.boundaryMesh()[patchi].start();
                    blockScalar* blk =
                        A_.interfaceCoeffs(blockIndex_[patchi]).data()
                      + pf*blockSize;
                    for (label k = 0; k < blockDim; ++k)
                    {
                        blk[blockP*blockDim + k] = 0;
                    }
                }
            }
        }

        // Row residual d (pRefValue - p_ref); the diagonal loop below adds
        // d p_ref to A x
        Ax_[c*blockDim + blockP] = 0;
        b_[c*blockDim + blockP] = d*pRefValue_;
    }

    // --- Row 3 diagonal contribution to A x
    if constexpr (zeroCopy)
    {
        // Zero-copy (D2.3): the diagonal blocks are the matrix array. Same
        // pass: row 3 scaled by s_p (in double, one rounding) and the clamp
        // count of all 16 entries (9.2 item 10)
        for (label celli = 0; celli < nCells; ++celli)
        {
            const doubleScalar x[blockDim] =
                {Ui[celli][0], Ui[celli][1], Ui[celli][2], pi[celli]};
            doubleScalar s = 0;
            for (label c = 0; c < blockDim; ++c)
            {
                s += Dd[di(celli, blockP, c)]*x[c];
            }
            Ax_[celli*blockDim + blockP] += s;

            diagScalar* __restrict__ blk = Dd + celli*blockSize;
            for (label k = 0; k < blockP*blockDim; ++k)
            {
                blk[k] = diagScalar(store(blk[k]));
            }
            for (label k = blockP*blockDim; k < blockSize; ++k)
            {
                blk[k] = diagScalar(store(sp*blk[k]));
            }
        }
    }
    else
    {
        for (label celli = 0; celli < nCells; ++celli)
        {
            const doubleScalar x[blockDim] =
                {Ui[celli][0], Ui[celli][1], Ui[celli][2], pi[celli]};
            doubleScalar s = 0;
            for (label c = 0; c < blockDim; ++c)
            {
                s += Dd[di(celli, blockP, c)]*x[c];
            }
            Ax_[celli*blockDim + blockP] += s;
        }

        // --- Narrow the diagonal blocks (row 3 scaled by s_p in double)
        blockScalarList& Ad = A_.diag();
        forAll(DdStage_, i)
        {
            const bool row3 = ((i % blockSize)/blockDim == blockP);
            Ad[i] = store(row3 ? sp*DdStage_[i] : DdStage_[i]);
        }
    }

    // --- Residual (right-hand side of the increment solve) and norms
    reduceScalar sums[6] = {0, 0, 0, 0, 0, 0};
    // 0: sum r^2 (all), 1: sum(|Ax|+|b|) all, 2: sum|r| rows 0-2,
    // 3: sum(|Ax|+|b|) rows 0-2, 4: sum|r| row 3, 5: sum(|Ax|+|b|) row 3
    const label refCell = (needRef_ ? pRefCell_ : -1);

    for (label celli = 0; celli < nCells; ++celli)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            const label i = celli*blockDim + k;
            const reduceScalar s = (k == blockP ? sp_ : 1.0);
            const reduceScalar r = s*(b_[i] - Ax_[i]);
            // The replaced continuity row of the pressure reference cell
            // (D-021) is not part of the physical system: it enters the
            // residual but not the normalisation sums
            const reduceScalar m =
                (k == blockP && celli == refCell)
              ? 0
              : s*(std::abs(Ax_[i]) + std::abs(b_[i]));

            rhs_[i] = store(r);

            sums[0] += r*r;
            sums[1] += m;
            if (k < blockP)
            {
                sums[2] += std::abs(r);
                sums[3] += m;
            }
            else
            {
                sums[4] += std::abs(r);
                sums[5] += m;
            }
        }
    }

    doubleReduce::parSum(sums, 6, mesh_.comm());

    // GUARD: normFactor >= cfVSmall (9.2, D3)
    normFactor_ = sums[1] + cfVSmall<reduceScalar>();
    residualL2_ = std::sqrt(max(sums[0], 0.0))/normFactor_;
    rU_ = sums[2]/(sums[3] + cfVSmall<reduceScalar>());
    rp_ = sums[4]/(sums[5] + cfVSmall<reduceScalar>());

    A_.markUpdated();

    if (timing_)
    {
        const doubleScalar tEnd = diagnostics::clock();
        times_.rhieChow += tr1 - tr0;
        times_.boundary += tb1 - tb0;
        times_.continuity += (tEnd - tStart) - (tr1 - tr0) - (tb1 - tb0);
    }
}


// ************************************************************************* //
