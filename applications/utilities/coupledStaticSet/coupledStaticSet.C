/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Application
    coupledStaticSet

Description
    Evaluates the coupledFoam static remediation set (spec 8.1 quality
    criteria + amendment C1 topological criteria) on the mesh of a case and
    prints the same once-per-run lines as coupledFoam. Read-only: nothing is
    written, so it can run on a cached mesh directory (e.g. run/T4a_mesh_np1)
    serial or decomposed.

    By default the library defaults (coupledDefaults.H) are used; with
    -fromCase the coupled.remediation entries of system/fvSolution.

Usage
    coupledStaticSet [-fromCase] [-parallel]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "remediation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Print the coupledFoam static remediation set counts (8.1 + C1),"
        " read-only"
    );
    argList::addBoolOption
    (
        "fromCase",
        "use coupled.remediation of system/fvSolution instead of the"
        " library defaults"
    );
    argList::noFunctionObjects();

    #include "setRootCase.H"
    #include "createTime.H"

    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    dictionary coupledDict;
    if (args.found("fromCase"))
    {
        coupledDict = mesh.solutionDict().subOrEmptyDict("coupled");
    }

    Info<< "coupledStaticSet: cells "
        << returnReduce(mesh.nCells(), sumOp<label>()) << ", ranks "
        << UPstream::nProcs() << ", solution dimensions "
        << mesh.nSolutionD() << nl << endl;

    remediation rem(mesh, coupledDict);
    rem.buildStatic();

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
