/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-doubleReduce

Description
    Documents why every reduction is accumulated in double (spec 6, 6.4).

    Sums N = 1e8 blockScalar (float) values of 1e-4, distributed over all
    ranks, once with a float accumulator and once with doubleReduce::sum.

    Pass (DECISIONS.md D-004): the double-accumulated sum equals the exact
    sum of the stored floats, N * double(float(1e-4)), to 1e-9 relative, and
    the float-accumulated sum is visibly wrong (> 1 % relative error).
    The distance of the exact sum to 1e4 (representation error of 1e-4 in
    float, 2.5e-8 relative) is reported for completeness.

    Field data (amendment D4, D-064): the same values stored as scalar
    (float in SP, double in DP), at most 1e7 per rank, through the
    templated doubleReduce::sum / average / weightedSum (weights 1): each
    must equal the exact sum (count) to 1e-9 relative.

Usage
    Test-doubleReduce [-parallel] [-json <file>] [-n <N>]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "doubleReduce.H"
#include "jsonWriter.H"
#include "PstreamReduceOps.H"
#include "scalarList.H"
#include <cmath>

using namespace Foam;

// Test parameters (fixed by spec 6.4)
static constexpr doubleScalar value = 1e-4;
static constexpr doubleScalar defaultN = 1e8;
static constexpr doubleScalar tolDouble = 1e-9;
static constexpr doubleScalar minFloatError = 1e-2;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("json", "file", "Write results as JSON");
    argList::addOption("n", "N", "Global number of values (1e8)");
    argList args(argc, argv);

    const doubleScalar nGlobalD = args.getOrDefault<doubleScalar>("n", defaultN);
    const label nGlobal = label(nGlobalD);
    const label nProcs = UPstream::nProcs();
    const label me = UPstream::myProcNo();

    label nLocal = nGlobal/nProcs;
    if (me == 0)
    {
        nLocal += nGlobal - nProcs*(nGlobal/nProcs);
    }

    const blockScalar v = narrow(value);
    blockScalarList a(nLocal, v);

    // Float accumulation (what the rule forbids)
    blockScalar floatLocal = 0;
    for (label i = 0; i < nLocal; ++i)
    {
        floatLocal += a[i];
    }
    // The cross-rank sum of the float partial sums is also done in float
    List<blockScalar> partials(nProcs, Zero);
    partials[me] = floatLocal;
    Foam::reduce
    (
        partials.data(),
        int(partials.size()),
        sumOp<blockScalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );
    blockScalar floatSum = 0;
    for (const blockScalar p : partials)
    {
        floatSum += p;
    }

    // Double accumulation
    const reduceScalar doubleSum = doubleReduce::sum(a);

    const reduceScalar exact = reduceScalar(nGlobal)*toDouble(v);
    const reduceScalar nominal = reduceScalar(nGlobal)*value;

    const reduceScalar relDouble = std::abs(doubleSum - exact)/exact;
    const reduceScalar relFloat = std::abs(toDouble(floatSum) - exact)/exact;
    const reduceScalar relRepr = std::abs(exact - nominal)/nominal;

    // Field data (D4): scalar lists through the templated reductions
    const label nField = min(nLocal, label(10000000));
    const scalarList af(nField, scalar(v));
    const scalarList wf(nField, scalar(1));
    const reduceScalar nFieldGlobal =
        doubleReduce::parSum(reduceScalar(nField), UPstream::worldComm);
    const reduceScalar exactField = nFieldGlobal*toDouble(v);
    const reduceScalar fieldSum = doubleReduce::sum(af);
    const reduceScalar fieldAvg = doubleReduce::average(af);
    const reduceScalar fieldWSum = doubleReduce::weightedSum(wf, af);
    const reduceScalar relField = max
    (
        std::abs(fieldSum - exactField)/exactField,
        max
        (
            std::abs(fieldAvg - toDouble(v))/toDouble(v),
            std::abs(fieldWSum - exactField)/exactField
        )
    );

    const bool pass =
        (relDouble < tolDouble) && (relFloat > minFloatError)
     && (relField < tolDouble);

    Info<< "N = " << nGlobal << " on " << nProcs << " ranks" << nl
        << "exact sum of stored floats = " << exact << nl
        << "double-accumulated sum     = " << doubleSum
        << "  rel. error " << relDouble << nl
        << "float-accumulated sum      = " << floatSum
        << "  rel. error " << relFloat << nl
        << "representation error vs " << nominal << " = " << relRepr << nl
        << "field data (D4): sum " << fieldSum << ", average " << fieldAvg
        << ", weighted sum " << fieldWSum << "  max rel. error " << relField
        << nl
        << (pass ? "PASS" : "FAIL") << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-doubleReduce");
        j.add("nProcs", nProcs);
        j.add("N", nGlobal);
        j.add("exactSumOfFloats", exact);
        j.add("doubleSum", doubleSum);
        j.add("floatSum", floatSum);
        j.add("relErrorDouble", relDouble);
        j.add("relErrorFloat", relFloat);
        j.add("relRepresentationError", relRepr);
        j.add("toleranceDouble", tolDouble);
        j.add("fieldN", label(nFieldGlobal));
        j.add("fieldSum", fieldSum);
        j.add("fieldAverage", fieldAvg);
        j.add("fieldWeightedSum", fieldWSum);
        j.add("relErrorField", relField);
        j.add("pass", pass);
        j.write(args.get<fileName>("json"));
    }

    return pass ? 0 : 1;
}


// ************************************************************************* //
