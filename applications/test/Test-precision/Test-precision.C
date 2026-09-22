/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-precision

Description
    Precision check per build (spec 3, 16; amendment D1; DECISIONS.md
    D-001, D-064).

    Reports the sizes of OpenFOAM's scalar and solveScalar and of the block
    library's blockScalar and reduceScalar and asserts, for the precision
    option the code was compiled with:

        build   sizeof(scalar)  sizeof(solveScalar)
        DP      8               8
        SP      4               4
        SPDP    4               8   (OpenFOAM's SPDP: fields float, linear
                                     solve double; the amendment's "8/4" is
                                     not an OpenFOAM build, D-001)

    and, in every build, sizeof(blockScalar) == 4 (the block matrix is
    float by design, D-001/D-064) and sizeof(reduceScalar) == 8.

Usage
    Test-precision [-json <file>]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "blockScalar.H"
#include "jsonWriter.H"
#include "foamVersion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noBanner();
    argList::addOption("json", "file", "Write results as JSON");
    argList args(argc, argv);

    const label sScalar = label(sizeof(scalar));
    const label sSolve = label(sizeof(solveScalar));
    const label sBlock = label(sizeof(blockScalar));
    const label sReduce = label(sizeof(reduceScalar));

    // Expected sizes of the build this binary was compiled for
#if defined(WM_SP)
    const char* build = "SP";
    const label expScalar = 4, expSolve = 4;
#elif defined(WM_SPDP)
    const char* build = "SPDP";
    const label expScalar = 4, expSolve = 8;
#elif defined(WM_DP)
    const char* build = "DP";
    const label expScalar = 8, expSolve = 8;
#else
    const char* build = "unknown";
    const label expScalar = -1, expSolve = -1;
#endif

    const bool buildOk = (sScalar == expScalar && sSolve == expSolve);
    const bool blockOk = (sBlock == 4 && sReduce == 8);
    const bool gate = (buildOk && blockOk);
    const bool literal = (sSolve == 4 && sScalar == 8);

    const char* opts = std::getenv("WM_OPTIONS");

    Info<< "WM_OPTIONS            = " << (opts ? opts : "unset") << nl
        << "precision option      = " << build << nl
        << "sizeof(scalar)        = " << sScalar
        << " (expected " << expScalar << ")" << nl
        << "sizeof(solveScalar)   = " << sSolve
        << " (expected " << expSolve << ")" << nl
        << "sizeof(blockScalar)   = " << sBlock << " (expected 4)" << nl
        << "sizeof(reduceScalar)  = " << sReduce << " (expected 8)" << nl
        << "gate (D1, D-064)      = " << (gate ? "PASS" : "FAIL") << nl
        << "spec literal criterion= " << (literal ? "true" : "false")
        << " (unsatisfiable, D-001)" << endl;

    if (args.found("json"))
    {
        jsonWriter j;
        j.add("test", "Test-precision");
        j.add("WM_OPTIONS", std::string(opts ? opts : "unset"));
        j.add("sizeof_scalar", sScalar);
        j.add("sizeof_solveScalar", sSolve);
        j.add("sizeof_blockScalar", sBlock);
        j.add("sizeof_reduceScalar", sReduce);
        j.add("precision", std::string(build));
        j.add("expected_sizeof_scalar", expScalar);
        j.add("expected_sizeof_solveScalar", expSolve);
        j.add("spec_literal_criterion", literal);
        j.add("pass", gate);
        j.write(args.get<fileName>("json"));
    }

    return gate ? 0 : 1;
}


// ************************************************************************* //
