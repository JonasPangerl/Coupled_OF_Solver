/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-nonOrthLimiter

Description
    Regression test of D-078: the limited non-orthogonal correction of
    nonOrthCorrection::correctionFromGrad.

    The limiter l = min(lambda |sn| / ((1 - lambda) |c|), 1) used to be
    evaluated as a quotient with the denominator guard cfVSmall = 1e-300.
    On a face with c = 0 exactly the quotient is lambda |sn| 1e300, which
    overflows once |sn| > ~1.8e8/lambda and raises SIGFPE under
    FOAM_SIGFPE (the F1 half-car, 2026-09-30). Two checks on the T0 cavity
    mesh, which keeps the correction path active because its round-off
    non-orthogonality (7.1e-14) exceeds orthogonalityTolerance 0:

    1. steep field p = 1e12 x: the correction is computed without a
       floating-point exception and is finite; the test also counts the
       faces on which the OLD formula would have overflowed and requires
       at least one, so that it really exercises the bug
    2. moderate field p = 1e-15 x + y (so that the limiter's quotient
       branch is taken on faces with k != 0): the result is bit-identical
       to the old formula on every face, as D-078 claims

Usage
    Test-nonOrthLimiter [-json <file>] [-parallel]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nonOrthCorrection.H"
#include "coupledConstants.H"
#include "linear.H"
#include "snGradScheme.H"
#include "zeroGradientFvPatchFields.H"
#include "jsonWriter.H"
#include <cmath>
#include <limits>

using namespace Foam;

// Test parameters (not solver parameters)
static constexpr scalar lambda = 0.5;          // nonOrthLimiter of the F1 case
static constexpr scalar steepSlope = 1e12;     // |sn| far above 1.8e8/lambda
static constexpr scalar tinySlope = 1e-15;     // sn << c on the x-faces

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static tmp<volScalarField> makeField
(
    const fvMesh& mesh,
    const word& name,
    const scalar ax,
    const scalar ay
)
{
    auto tp = tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField& p = tp.ref();
    const volVectorField& C = mesh.C();
    forAll(p, celli)
    {
        p[celli] = ax*C[celli].x() + ay*C[celli].y();
    }
    p.correctBoundaryConditions();
    return tp;
}


// The formula before D-078, internal faces only (the ones that matter:
// boundary faces of the cavity have k = 0 and sn = 0 with zeroGradient).
// Returns the corrected field; overflow faces are only COUNTED here, never
// evaluated, so the test itself cannot trap.
static scalarField oldFormula
(
    const volScalarField& vf,
    const volVectorField& gradVf,
    label& nOverflow
)
{
    const fvMesh& mesh = vf.mesh();
    const tmp<surfaceScalarField> tc =
        linear<vector>(mesh).dotInterpolate
        (
            mesh.nonOrthCorrectionVectors(),
            gradVf
        );
    const tmp<surfaceScalarField> tsn =
        fv::snGradScheme<scalar>::snGrad
        (
            vf,
            tmp<surfaceScalarField>(mesh.nonOrthDeltaCoeffs()),
            "SndGrad"
        );
    const scalarField& c = tc().primitiveField();
    const scalarField& sn = tsn().primitiveField();

    const scalar big = std::numeric_limits<scalar>::max();
    scalarField out(c.size());
    nOverflow = 0;
    forAll(c, facei)
    {
        const scalar num = lambda*mag(sn[facei]);
        const scalar den = (1 - lambda)*mag(c[facei]) + cfVSmall<scalar>();
        // num/den > big  <=>  num > big*den (big*den may itself be inf,
        // then the comparison is false, correctly)
        if (num > big*den)
        {
            ++nOverflow;
            out[facei] = c[facei];     // min(inf, 1) = 1
        }
        else
        {
            out[facei] = c[facei]*min(num/den, scalar(1));
        }
    }
    return out;
}


int main(int argc, char *argv[])
{
    argList::addOption("json", "file", "Write results as JSON");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // orthogonalityTolerance 0: the exact test, as the solver default
    const nonOrthCorrection noc(mesh, lambda, lambda, 0);

    const label nFaces = returnReduce(mesh.nInternalFaces(), sumOp<label>());

    // --- 1. steep gradient: no SIGFPE (a trap would abort the process)
    label nOverflowOld = 0;
    label nNonFinite = 0;
    {
        const tmp<volScalarField> tp = makeField(mesh, "pSteep", steepSlope, 0);
        const volVectorField g(fvc::grad(tp()));
        const tmp<surfaceScalarField> tcorr = noc.correctionFromGrad(tp(), g);
        for (const scalar v : tcorr().primitiveField())
        {
            nNonFinite += (std::isfinite(v) ? 0 : 1);
        }
        oldFormula(tp(), g, nOverflowOld);
    }
    reduce(nOverflowOld, sumOp<label>());
    reduce(nNonFinite, sumOp<label>());

    // --- 2. moderate field: bit-identical to the old formula
    label nMismatch = 0;
    label nQuotient = 0;
    {
        const tmp<volScalarField> tp =
            makeField(mesh, "pModerate", tinySlope, 1);
        const volVectorField g(fvc::grad(tp()));
        const tmp<surfaceScalarField> tcorr = noc.correctionFromGrad(tp(), g);
        label nOv = 0;
        const scalarField ref(oldFormula(tp(), g, nOv));
        const scalarField& c = tcorr().primitiveField();
        forAll(c, facei)
        {
            // Bitwise: == on doubles, and both NaN counts as equal
            const bool same =
                (c[facei] == ref[facei])
             || (std::isnan(c[facei]) && std::isnan(ref[facei]));
            nMismatch += (same ? 0 : 1);
        }
        // Faces where the quotient branch was taken (the limiter reduced c)
        const tmp<surfaceScalarField> tc =
            linear<vector>(mesh).dotInterpolate
            (
                mesh.nonOrthCorrectionVectors(),
                g
            );
        const scalarField& cRaw = tc().primitiveField();
        forAll(c, facei)
        {
            nQuotient += (mag(c[facei]) < mag(cRaw[facei]) ? 1 : 0);
        }
        nMismatch += nOv;      // the moderate field must not overflow at all
    }
    reduce(nMismatch, sumOp<label>());
    reduce(nQuotient, sumOp<label>());

    const bool pass =
        !noc.orthogonal()
     && nOverflowOld > 0
     && nNonFinite == 0
     && nMismatch == 0
     && nQuotient > 0;

    Info<< "internal faces " << nFaces << "  ranks " << UPstream::nProcs()
        << nl
        << "correction path active (not orthogonal): " << !noc.orthogonal()
        << nl
        << "1. steep p = " << steepSlope << " x: faces the old formula"
        << " overflows on " << nOverflowOld << ", non-finite results "
        << nNonFinite << nl
        << "2. moderate p: faces limited by the quotient " << nQuotient
        << ", bitwise mismatches vs old formula " << nMismatch << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-nonOrthLimiter");
        j.add("nProcs", UPstream::nProcs());
        j.add("nInternalFaces", nFaces);
        j.add("correctionActive", !noc.orthogonal());
        j.add("nOverflowOld", nOverflowOld);
        j.add("nNonFinite", nNonFinite);
        j.add("nQuotient", nQuotient);
        j.add("nMismatch", nMismatch);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
