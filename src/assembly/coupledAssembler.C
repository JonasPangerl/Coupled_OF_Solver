/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledAssembler.H"
#include "coupledDefaults.H"
#include "doubleReduce.H"
#include "fvCFD.H"
#include "laplacianScheme.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "wedgeFvPatch.H"
#include "findRefCell.H"
#include "IStringStream.H"

using namespace Foam::boundaryCoupling;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::coupledAssembler::zeroFluxPatch(const fvPatch& patch)
{
    return
        isA<symmetryFvPatch>(patch)
     || isA<symmetryPlaneFvPatch>(patch)
     || isA<wedgeFvPatch>(patch);
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
    rc_(mesh),
    mrfPtr_(mrfPtr),
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
    A_(mesh),
    rhs_(blockDim*mesh.nCells(), Zero),
    Dd_(blockSize*mesh.nCells(), Zero),
    Ax_(blockDim*mesh.nCells(), Zero),
    b_(blockDim*mesh.nCells(), Zero),
    rMom_(mesh.nCells(), Zero),
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledAssembler::setStaticCells(const boolList& isStatic)
{
    noc_.setStaticCells(isStatic);
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
    UPtr_ = &U;
    pPtr_ = &p;

    A_.clear();
    Dd_ = Zero;
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
            Dd_[di(celli, c, c)] += diagU[celli];
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
                    Dd_[di(celli, r, c)] += C[celli][r*blockP + c];
                }
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
            Dd_[di(P, c, blockP)] += wf*S[c];
            const doubleScalar ug = (1 - wf)*S[c];
            Ub[c*blockDim + blockP] = store(ug);
            Ax_[P*blockDim + c] += ug*pi[N];

            Dd_[di(N, c, blockP)] -= (1 - wf)*S[c];
            const doubleScalar lg = -wf*S[c];
            Lb[c*blockDim + blockP] = store(lg);
            Ax_[N*blockDim + c] += lg*pi[P];
        }
    }

    // --- Boundary faces
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
                        Dd_[di(P, c, c)] += ic[pf][c];
                        // Native sign: result -= bouCoeffs*psi_nbr
                        const doubleScalar cc = -bc[pf][c];
                        Cb[c*blockDim + c] = store(cc);
                        Ax_[P*blockDim + c] += cc*Un[pf][c];

                        Dd_[di(P, c, blockP)] += wp[pf]*Sp[pf][c];
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
                        Dd_[di(P, c, c)] += ic[pf][c];
                        b_[P*blockDim + c] += bc[pf][c]*Un[pf][c];

                        Dd_[di(P, c, blockP)] += wp[pf]*Sp[pf][c];
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
                        Dd_[di(P, c, c)] += ic[pf][c];
                        b_[P*blockDim + c] += bc[pf][c];

                        // p_f = viC p_P + vbC
                        Dd_[di(P, c, blockP)] += Sp[pf][c]*viCp[pf];
                        b_[P*blockDim + c] -= Sp[pf][c]*vbCp[pf];
                    }
                }
                break;
            }
        }
    }

    // --- Diagonal contribution to A x for rows 0-2 (PTC excluded: it
    //     cancels in b - A x), momentum residual
    forAll(V, celli)
    {
        const doubleScalar x[blockDim] =
            {Ui[celli][0], Ui[celli][1], Ui[celli][2], pi[celli]};

        for (label r = 0; r < blockP; ++r)
        {
            doubleScalar s = 0;
            for (label c = 0; c < blockDim; ++c)
            {
                s += Dd_[di(celli, r, c)]*x[c];
            }
            Ax_[celli*blockDim + r] += s;
            rMom_[celli][r] = b_[celli*blockDim + r] - Ax_[celli*blockDim + r];
        }
    }
}


void Foam::coupledAssembler::assembleContinuity(const scalarField& rDeltaTV)
{
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

    // --- PTC diagonal and abar
    for (label celli = 0; celli < nCells; ++celli)
    {
        doubleScalar s = 0;
        for (label c = 0; c < blockP; ++c)
        {
            Dd_[di(celli, c, c)] += rDeltaTV[celli];
            s += Dd_[di(celli, c, c)];
        }
        abar_[celli] = s/blockP;
    }

    // --- Rhie-Chow D, D_f, explicit face term q
    rc_.updateD(abar_, p);
    rc_.updateExplicit(p, noc_);
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
            Dd_[di(P, blockP, c)] += wf*S[c];
            const doubleScalar ud = (1 - wf)*S[c];
            Ub[blockP*blockDim + c] = store(ud);
            Ax_[P*blockDim + blockP] += ud*Ui[N][c];

            // Neighbour row (-)
            Dd_[di(N, blockP, c)] -= (1 - wf)*S[c];
            const doubleScalar ld = -wf*S[c];
            Lb[blockP*blockDim + c] = store(ld);
            Ax_[N*blockDim + blockP] += ld*Ui[P][c];
        }

        Dd_[di(P, blockP, blockP)] += g;
        Ub[blockP*blockDim + blockP] = store(-g);
        Ax_[P*blockDim + blockP] -= g*pi[N];

        Dd_[di(N, blockP, blockP)] += g;
        Lb[blockP*blockDim + blockP] = store(-g);
        Ax_[N*blockDim + blockP] -= g*pi[P];

        b_[P*blockDim + blockP] -= q[facei];
        b_[N*blockDim + blockP] += q[facei];
    }

    // --- Boundary faces, row 3
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
                        Dd_[di(P, blockP, c)] += wp[pf]*Sp[pf][c];
                        const doubleScalar cd = (1 - wp[pf])*Sp[pf][c];
                        Cb[blockP*blockDim + c] = store(cd);
                        Ax_[P*blockDim + blockP] += cd*Un[pf][c];
                    }

                    Dd_[di(P, blockP, blockP)] += g;
                    Cb[blockP*blockDim + blockP] = store(-g);
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
                        Dd_[di(P, blockP, c)] += wp[pf]*Sp[pf][c];
                    }
                    b_[P*blockDim + blockP] -=
                        (1 - wp[pf])*(Sp[pf] & Un[pf]);

                    Dd_[di(P, blockP, blockP)] += g;
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
                        Dd_[di(P, blockP, c)] += Sp[pf][c]*viCU[pf][c];
                    }
                    b_[P*blockDim + blockP] -= Sp[pf] & vbCU[pf];

                    // -g (p_b - p_P), p_b = viC p_P + vbC (g = 0 unless
                    // p is fixed on this patch, D-013)
                    const doubleScalar g = rc_.g(patchi, pf);
                    Dd_[di(P, blockP, blockP)] += g*(1 - viCp[pf]);
                    b_[P*blockDim + blockP] += g*vbCp[pf];

                    b_[P*blockDim + blockP] -= qp[pf];
                }
                break;
            }
        }
    }

    // --- Pressure reference (native setReference: diag doubled, the
    //     residual drives p_ref towards pRefValue)
    if (needRef_ && pRefCell_ >= 0)
    {
        const doubleScalar d = Dd_[di(pRefCell_, blockP, blockP)];
        Dd_[di(pRefCell_, blockP, blockP)] += d;
        Ax_[pRefCell_*blockDim + blockP] += d*pi[pRefCell_];
        b_[pRefCell_*blockDim + blockP] += d*pRefValue_;
    }

    // --- Row 3 diagonal contribution to A x
    for (label celli = 0; celli < nCells; ++celli)
    {
        const doubleScalar x[blockDim] =
            {Ui[celli][0], Ui[celli][1], Ui[celli][2], pi[celli]};
        doubleScalar s = 0;
        for (label c = 0; c < blockDim; ++c)
        {
            s += Dd_[di(celli, blockP, c)]*x[c];
        }
        Ax_[celli*blockDim + blockP] += s;
    }

    // --- Narrow the diagonal blocks
    blockScalarList& Ad = A_.diag();
    forAll(Dd_, i)
    {
        Ad[i] = store(Dd_[i]);
    }

    // --- Continuity row scaling (5.3f); GUARD: s_p <= 1/VSMALL
    sp_ = 1.0/max(rc_.Dref(), VSMALL);
    A_.scaleRow(blockP, narrow(sp_));

    // --- Residual (right-hand side of the increment solve) and norms
    reduceScalar sums[6] = {0, 0, 0, 0, 0, 0};
    // 0: sum r^2 (all), 1: sum(|Ax|+|b|) all, 2: sum|r| rows 0-2,
    // 3: sum(|Ax|+|b|) rows 0-2, 4: sum|r| row 3, 5: sum(|Ax|+|b|) row 3

    for (label celli = 0; celli < nCells; ++celli)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            const label i = celli*blockDim + k;
            const reduceScalar s = (k == blockP ? sp_ : 1.0);
            const reduceScalar r = s*(b_[i] - Ax_[i]);
            const reduceScalar m = s*(std::abs(Ax_[i]) + std::abs(b_[i]));

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

    doubleReduce::parSum(sums, 6, UPstream::worldComm);

    // GUARD: normFactor >= SMALL (9.2)
    normFactor_ = sums[1] + SMALL;
    residualL2_ = std::sqrt(max(sums[0], 0.0))/normFactor_;
    rU_ = sums[2]/(sums[3] + SMALL);
    rp_ = sums[4]/(sums[5] + SMALL);

    A_.markUpdated();
}


// ************************************************************************* //
