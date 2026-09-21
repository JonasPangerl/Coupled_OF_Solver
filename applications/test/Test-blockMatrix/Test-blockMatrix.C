/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-blockMatrix

Description
    Unit test of blockLduMatrix4::Amul (spec 6.4).

    Builds the 4x4-block Poisson-like system of testPoissonSystem.H on the
    case mesh and compares the block Amul, component by component, with the
    native lduMatrix::Amul of the scalar operator (including processor
    interfaces in parallel).

    Pass: max_k,i |(A x)_block - (A x)_native| < 1e-6 * ||x||_2 (global).
    Written for the 1-vs-N-rank comparison (done by the pytest suite, 1e-6
    relative): per-component global sums of Ax and |Ax|.

Usage
    Test-blockMatrix [-parallel] [-json <file>]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "testPoissonSystem.H"
#include "doubleReduce.H"
#include "jsonWriter.H"

using namespace Foam;

// Test tolerance (fixed by spec 6.4)
static constexpr doubleScalar tolAmul = 1e-6;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("json", "file", "Write results as JSON");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    tmp<volScalarField> tpsi = testPoisson::makePsi(mesh);
    volScalarField& psi = tpsi.ref();
    tmp<fvScalarMatrix> tM = testPoisson::makeScalarMatrix(psi);
    const fvScalarMatrix& M = tM();

    blockLduMatrix4 B(mesh);
    testPoisson::fillBlockMatrix(M, B);

    const blockScalarList x = testPoisson::makeX(mesh);
    blockScalarList Ax(B.nRows());
    B.Amul(Ax, x);

    // Native operator: lduMatrix with the boundary-including diagonal
    lduMatrix nat(M);
    nat.diag() = M.D();
    const lduInterfaceFieldPtrsList interfaces =
        psi.boundaryField().scalarInterfaces();

    const label nCells = mesh.nCells();
    doubleScalar maxDiff = 0;
    FixedList<doubleScalar, 4> sumAx(Zero);
    FixedList<doubleScalar, 4> sumMagAx(Zero);

    for (label k = 0; k < blockDim; ++k)
    {
        solveScalarField xk(nCells);
        for (label celli = 0; celli < nCells; ++celli)
        {
            xk[celli] = solveScalar(x[celli*blockDim + k]);
        }

        solveScalarField Axk(nCells, Zero);
        nat.Amul
        (
            Axk,
            tmp<solveScalarField>(new solveScalarField(xk)),
            M.boundaryCoeffs(),
            interfaces,
            0
        );

        for (label celli = 0; celli < nCells; ++celli)
        {
            const doubleScalar b = toDouble(Ax[celli*blockDim + k]);
            maxDiff = std::max
            (
                maxDiff,
                std::abs(b - doubleScalar(Axk[celli]))
            );
            sumAx[k] += b;
            sumMagAx[k] += std::abs(b);
        }
    }

    Foam::reduce(maxDiff, maxOp<doubleScalar>());
    doubleReduce::parSum(sumAx.data(), 4, UPstream::worldComm);
    doubleReduce::parSum(sumMagAx.data(), 4, UPstream::worldComm);

    const doubleScalar xNorm = doubleReduce::norm2(x);
    const bool pass = maxDiff < tolAmul*xNorm;

    Info<< "cells " << returnReduce(nCells, sumOp<label>())
        << "  ranks " << UPstream::nProcs() << nl
        << "||x||_2 = " << xNorm << nl
        << "max |Ax_block - Ax_native| = " << maxDiff
        << "  (limit " << tolAmul*xNorm << ")" << nl
        << "sum(Ax)   per component " << sumAx << nl
        << "sum|Ax|   per component " << sumMagAx << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-blockMatrix");
        j.add("nProcs", UPstream::nProcs());
        j.add("nCells", returnReduce(nCells, sumOp<label>()));
        j.add("xNorm2", xNorm);
        j.add("maxAbsDiff", maxDiff);
        j.add("limit", tolAmul*xNorm);
        j.addList("sumAx", List<doubleScalar>(sumAx));
        j.addList("sumMagAx", List<doubleScalar>(sumMagAx));
        j.add("nnzBlocks", returnReduce(B.nnzBlocks(), sumOp<label>()));
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
