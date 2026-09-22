/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "precisionProfile.H"
#include "coupledDefaults.H"
#include "error.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
    Foam::precisionProfile::values dpValues()
    {
        namespace d = Foam::coupledDefaults::dpProfile;
        return
        {
            "dp", d::residualTol, d::tolerance, d::etaMin, d::boundOmegaMin,
            d::kMin, d::andersonAboveMaxCells, d::andersonLargeMaxM
        };
    }

    Foam::precisionProfile::values spValues()
    {
        namespace d = Foam::coupledDefaults::spProfile;
        return
        {
            "sp", d::residualTol, d::tolerance, d::etaMin, d::boundOmegaMin,
            d::kMin, d::andersonAboveMaxCells, d::andersonLargeMaxM
        };
    }
}


// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

Foam::precisionProfile::values Foam::precisionProfile::current_ =
    (sizeof(Foam::scalar) == 4 ? spValues() : dpValues());

Foam::word Foam::precisionProfile::keyword_("auto");


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::precisionProfile::autoName()
{
    return (sizeof(scalar) == 4 ? word("sp") : word("dp"));
}


Foam::precisionProfile::values Foam::precisionProfile::named
(
    const word& name
)
{
    if (name == "dp")
    {
        return dpValues();
    }
    if (name == "sp")
    {
        return spValues();
    }
    FatalErrorInFunction
        << "precisionProfile must be auto, dp or sp, not " << name
        << exit(FatalError);
    return dpValues();
}


const Foam::precisionProfile::values& Foam::precisionProfile::select
(
    const dictionary& coupledDict
)
{
    keyword_ = coupledDict.getOrDefault<word>
    (
        "precisionProfile",
        word(coupledDefaults::precisionProfile)
    );
    if (keyword_ != "auto" && keyword_ != "dp" && keyword_ != "sp")
    {
        FatalIOErrorInFunction(coupledDict)
            << "coupled.precisionProfile must be auto, dp or sp, not "
            << keyword_ << exit(FatalIOError);
    }
    current_ = named(keyword_ == "auto" ? autoName() : keyword_);

    const values& v = current_;
    Info<< "coupledFoam: precision profile " << v.name << " ("
        << keyword_ << ", sizeof(scalar) " << label(sizeof(scalar))
        << "): defaults residualTol " << v.residualTol
        << ", tolerance " << v.tolerance
        << ", etaMin " << v.etaMin
        << ", bounds.omegaMin " << v.boundOmegaMin
        << ", bounds.kMin " << v.kMin
        << ", Anderson above anderson.maxCells "
        << (v.andersonAboveMaxCells ? "allowed (m <= " : "disabled")
        << (v.andersonAboveMaxCells ? Foam::name(v.andersonLargeMaxM) : "")
        << (v.andersonAboveMaxCells ? ")" : "")
        << "; explicitly set keywords win" << endl;

    return current_;
}


void Foam::precisionProfile::writeSettings(dictionary& dict)
{
    const values& v = current_;
    dictionary d;
    d.add("keyword", keyword_);
    d.add("profile", v.name);
    d.add("sizeofScalar", label(sizeof(scalar)));
    d.add("residualTol", v.residualTol);
    d.add("tolerance", v.tolerance);
    d.add("etaMin", v.etaMin);
    d.add("boundOmegaMin", v.boundOmegaMin);
    d.add("kMin", v.kMin);
    d.add("andersonAboveMaxCells", v.andersonAboveMaxCells);
    d.add("andersonLargeMaxM", v.andersonLargeMaxM);
    dict.add("precisionProfile", d);
}


// ************************************************************************* //
