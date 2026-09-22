/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledForcesFO.H"
#include "turbulentTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(coupledForcesFO, 0);
    addToRunTimeSelectionTable(functionObject, coupledForcesFO, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::coupledForcesFO::coupledForcesFO
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    forces_(),
    dormantReported_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::coupledForcesFO::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    if (dict.found("patches"))
    {
        forces_.reset(new coupledForces(mesh_, dict, name()));
    }
    else
    {
        forces_.reset
        (
            new coupledForces(mesh_, coupledForces::readDict(time_), name())
        );
    }
    return true;
}


bool Foam::functionObjects::coupledForcesFO::execute()
{
    if (coupledForces::internal())
    {
        if (!dormantReported_)
        {
            Info<< type() << ' ' << name() << ": dormant, coupledFoam"
                << " evaluates coupledForces after every outer iteration"
                << endl;
            dormantReported_ = true;
        }
        return true;
    }

    typedef incompressible::turbulenceModel icoModel;
    const auto* turbp = mesh_.cfindObject<icoModel>(icoModel::propertiesName);
    if (!turbp)
    {
        FatalErrorInFunction
            << type() << ' ' << name() << ": no incompressible turbulence"
            << " model (" << icoModel::propertiesName << ") in the database"
            << exit(FatalError);
    }

    forces_->evaluate(turbp->devReff()());

    const coupledForces::coefficients& c = forces_->coeffs();
    setResult("Cd", scalar(c.Cd));
    setResult("Cl", scalar(c.Cl));
    setResult("CmPitch", scalar(c.Cm));
    setResult("Cd_p", scalar(c.Cd_p));
    setResult("Cd_v", scalar(c.Cd_v));
    setResult("Cl_p", scalar(c.Cl_p));
    setResult("Cl_v", scalar(c.Cl_v));

    return true;
}


bool Foam::functionObjects::coupledForcesFO::write()
{
    if (!coupledForces::internal() && forces_ && forces_->evaluated())
    {
        forces_->write(double(time_.timeOutputValue()));
    }
    return true;
}


// ************************************************************************* //
