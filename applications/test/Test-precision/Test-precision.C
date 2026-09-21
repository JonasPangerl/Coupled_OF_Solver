/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    Test-precision

Description
    Phase 0 precision check (spec 3, 16; DECISIONS.md D-001).

    Reports the sizes of OpenFOAM's scalar and solveScalar and of the block
    library's blockScalar and reduceScalar. Gate criterion (D-001):
        sizeof(blockScalar) == 4, sizeof(reduceScalar) == 8,
        sizeof(scalar) == 8.
    The spec's literal criterion (sizeof(solveScalar) == 4 together with
    sizeof(scalar) == 8) is reported as well; it is unsatisfiable with any
    OpenFOAM precision option.

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

    const bool gate = (sBlock == 4 && sReduce == 8 && sScalar == 8);
    const bool literal = (sSolve == 4 && sScalar == 8);

    const char* opts = std::getenv("WM_OPTIONS");

    Info<< "WM_OPTIONS            = " << (opts ? opts : "unset") << nl
        << "sizeof(scalar)        = " << sScalar << nl
        << "sizeof(solveScalar)   = " << sSolve << nl
        << "sizeof(blockScalar)   = " << sBlock << nl
        << "sizeof(reduceScalar)  = " << sReduce << nl
        << "gate (D-001)          = " << (gate ? "PASS" : "FAIL") << nl
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
        j.add("spec_literal_criterion", literal);
        j.add("pass", gate);
        j.write(args.get<fileName>("json"));
    }

    return gate ? 0 : 1;
}


// ************************************************************************* //
