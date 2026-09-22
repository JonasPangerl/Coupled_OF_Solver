/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledForces.H"
#include "IOdictionary.H"
#include "surfaceFields.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "doubleVector.H"
#include <cmath>
#include <cstdio>
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::coupledForces::dictName("coupledForcesDict");

bool Foam::coupledForces::internal_ = false;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Guard of the double normalisations and reference scales
constexpr double tinyD = 1e-300;   // D3: switch to cfVSmall<double>() on merge

inline double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

inline Foam::vector toVector(const double a[3])
{
    return Foam::vector
    (
        Foam::scalar(a[0]),
        Foam::scalar(a[1]),
        Foam::scalar(a[2])
    );
}

inline void readVec
(
    const Foam::dictionary& dict,
    const Foam::word& key,
    double v[3]
)
{
    const Foam::doubleVector d(dict.get<Foam::doubleVector>(key));
    v[0] = d.x();
    v[1] = d.y();
    v[2] = d.z();
}

std::string fmtVec(const double v[3])
{
    char buf[128];
    std::snprintf(buf, sizeof(buf), "(%.10g %.10g %.10g)", v[0], v[1], v[2]);
    return buf;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coupledForces::readSettings(const dictionary& dict)
{
    patchIDs_ =
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches")).sortedToc();

    if (returnReduceOr(patchIDs_.empty()))
    {
        FatalIOErrorInFunction(dict)
            << "coupledForces: no patch matches " << dict.get<wordRes>("patches")
            << exit(FatalIOError);
    }

    pName_ = dict.getOrDefault<word>("p", "p");
    rhoInf_ = dict.get<doubleScalar>("rhoInf");
    magUInf_ = dict.get<doubleScalar>("magUInf");
    lRef_ = dict.get<doubleScalar>("lRef");
    Aref_ = dict.get<doubleScalar>("Aref");
    pRef_ = dict.getOrDefault<doubleScalar>("pRef", 0);

    if
    (
        !(rhoInf_ > tinyD) || !(magUInf_ > tinyD)
     || !(lRef_ > tinyD) || !(Aref_ > tinyD)
    )
    {
        FatalIOErrorInFunction(dict)
            << "coupledForces: rhoInf, magUInf, lRef and Aref must be > 0"
            << exit(FatalIOError);
    }

    readVec(dict, "CofR", CofR_);

    // Frame of native forceCoeffs: cartesian(CofR, e3 = liftDir,
    // e1 = dragDir), coordinateRotations::axes::rotation (E3_E1)
    double lift[3], drag[3];
    readVec(dict, "liftDir", lift);
    readVec(dict, "dragDir", drag);

    const double ml = std::sqrt(dot3(lift, lift));
    if (!(ml > tinyD))
    {
        FatalIOErrorInFunction(dict)
            << "coupledForces: liftDir has zero length" << exit(FatalIOError);
    }
    for (int i = 0; i < 3; ++i)
    {
        e3_[i] = lift[i]/ml;
    }
    const double de = dot3(drag, e3_);
    for (int i = 0; i < 3; ++i)
    {
        e1_[i] = drag[i] - de*e3_[i];
    }
    const double md = std::sqrt(dot3(e1_, e1_));
    if (!(md > 1e-6*std::sqrt(dot3(drag, drag))) || !(md > tinyD))
    {
        FatalIOErrorInFunction(dict)
            << "coupledForces: dragDir is zero or parallel to liftDir"
            << exit(FatalIOError);
    }
    for (int i = 0; i < 3; ++i)
    {
        e1_[i] /= md;
    }
    cross3(e3_, e1_, e2_);

    if (dict.found("pitchAxis"))
    {
        double pa[3];
        readVec(dict, "pitchAxis", pa);
        const double mp = std::sqrt(dot3(pa, pa));
        const double c = (mp > tinyD ? dot3(pa, e2_)/mp : 0);
        if (c < 1 - 1e-6)
        {
            Info<< "coupledForces " << name_ << ": pitchAxis "
                << fmtVec(pa).c_str() << " differs from the native pitch axis"
                << " e2 = liftDir ^ dragDir = " << fmtVec(e2_).c_str()
                << "; Cm is taken about e2 like native forceCoeffs (v2606"
                << " ignores pitchAxis)" << endl;
        }
    }
}


void Foam::coupledForces::openFile()
{
    if (!file_.empty() || !UPstream::master())
    {
        return;
    }

    const Time& t = mesh_.time();
    const fileName dir
    (
        t.globalPath()/"postProcessing"/name_
       /t.timeName(t.startTime().value())
    );
    mkDir(dir);
    file_ = dir/"coeffs.dat";

    const bool fresh = !isFile(file_) || fileSize(file_) <= 0;
    if (!fresh)
    {
        return;
    }

    std::ofstream os(file_.c_str(), std::ios::out | std::ios::app);
    char buf[256];
    os  << "# coupledForces (amendment D6): force coefficients, double"
        << " accumulation" << '\n';
    os  << "# patches";
    for (const label patchi : patchIDs_)
    {
        os  << ' ' << mesh_.boundary()[patchi].name();
    }
    os  << '\n';
    std::snprintf
    (
        buf, sizeof(buf),
        "# rhoInf %.10g magUInf %.10g lRef %.10g Aref %.10g pRef %.10g\n",
        rhoInf_, magUInf_, lRef_, Aref_, pRef_
    );
    os  << buf;
    os  << "# CofR " << fmtVec(CofR_) << " (mesh frame)" << '\n';
    os  << "# dragDir e1 " << fmtVec(e1_) << " liftDir e3 " << fmtVec(e3_)
        << " pitch axis e2 = e3 ^ e1 " << fmtVec(e2_) << '\n';
    os  << "# iter Cd Cl Cm Cd_p Cd_v Cl_p Cl_v" << '\n';
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledForces::coupledForces
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    mesh_(mesh),
    name_(name),
    patchIDs_(),
    pName_("p"),
    rhoInf_(1),
    magUInf_(1),
    lRef_(1),
    Aref_(1),
    pRef_(0),
    CofR_{0, 0, 0},
    e1_{1, 0, 0},
    e2_{0, 1, 0},
    e3_{0, 0, 1},
    Fp_{0, 0, 0},
    Fv_{0, 0, 0},
    Mp_{0, 0, 0},
    Mv_{0, 0, 0},
    c_(),
    evaluated_(false),
    file_()
{
    readSettings(dict);

    Info<< "coupledForces " << name_ << ": patches";
    for (const label patchi : patchIDs_)
    {
        Info<< ' ' << mesh_.boundary()[patchi].name();
    }
    Info<< ", CofR " << fmtVec(CofR_).c_str()
        << ", e1 (drag) " << fmtVec(e1_).c_str()
        << ", e3 (lift) " << fmtVec(e3_).c_str()
        << ", e2 (pitch) " << fmtVec(e2_).c_str() << endl;
}


Foam::autoPtr<Foam::coupledForces> Foam::coupledForces::New
(
    const fvMesh& mesh,
    const word& name
)
{
    if (!dictPresent(mesh.time()))
    {
        return nullptr;
    }
    return autoPtr<coupledForces>::New(mesh, readDict(mesh.time()), name);
}


bool Foam::coupledForces::dictPresent(const Time& runTime)
{
    IOobject io
    (
        dictName,
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
    return io.typeHeaderOk<IOdictionary>(true);
}


Foam::dictionary Foam::coupledForces::readDict(const Time& runTime)
{
    IOdictionary d
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    return dictionary(d);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledForces::evaluate
(
    const volScalarField& p,
    const volSymmTensorField& devReff
)
{
    // Kinematic pressure: the native incompressible scaling (rho rhoInf),
    // pRef given in dynamic units as in native forces
    const double rho = rhoInf_;
    const double pRefKin =
        (p.dimensions() == dimPressure ? pRef_ : pRef_/rhoInf_);
    const double rhoP = (p.dimensions() == dimPressure ? 1.0 : rhoInf_);

    // 0-2 Fp, 3-5 Fv, 6-8 Mp, 9-11 Mv
    double s[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    const auto& Sfb = mesh_.Sf().boundaryField();
    const auto& Cb = mesh_.C().boundaryField();
    const auto& pb = p.boundaryField();
    const auto& Rb = devReff.boundaryField();

    for (const label patchi : patchIDs_)
    {
        const vectorField& Sf = Sfb[patchi];
        const vectorField& Cf = Cb[patchi];
        const scalarField& pf = pb[patchi];
        const symmTensorField& R = Rb[patchi];

        forAll(Sf, facei)
        {
            const double S[3] =
            {
                double(Sf[facei].x()),
                double(Sf[facei].y()),
                double(Sf[facei].z())
            };
            const double d[3] =
            {
                double(Cf[facei].x()) - CofR_[0],
                double(Cf[facei].y()) - CofR_[1],
                double(Cf[facei].z()) - CofR_[2]
            };
            const symmTensor& r = R[facei];
            const double xx = double(r.xx()), xy = double(r.xy());
            const double xz = double(r.xz()), yy = double(r.yy());
            const double yz = double(r.yz()), zz = double(r.zz());

            const double pp = rhoP*(double(pf[facei]) - pRefKin);
            const double fp[3] = {pp*S[0], pp*S[1], pp*S[2]};

            // S . devReff (symmetric), times rho
            const double fv[3] =
            {
                rho*(S[0]*xx + S[1]*xy + S[2]*xz),
                rho*(S[0]*xy + S[1]*yy + S[2]*yz),
                rho*(S[0]*xz + S[1]*yz + S[2]*zz)
            };

            double mp[3], mv[3];
            cross3(d, fp, mp);
            cross3(d, fv, mv);

            for (int i = 0; i < 3; ++i)
            {
                s[i] += fp[i];
                s[3 + i] += fv[i];
                s[6 + i] += mp[i];
                s[9 + i] += mv[i];
            }
        }
    }

    // One MPI_DOUBLE reduction of the 12 sums (D2.4)
    if (UPstream::parRun())
    {
        Foam::reduce(s, 12, sumOp<double>(), UPstream::msgType());
    }

    for (int i = 0; i < 3; ++i)
    {
        Fp_[i] = s[i];
        Fv_[i] = s[3 + i];
        Mp_[i] = s[6 + i];
        Mv_[i] = s[9 + i];
    }

    const double q = 0.5*rhoInf_*magUInf_*magUInf_;
    const double fs = 1.0/std::max(q*Aref_, tinyD);
    const double ms = 1.0/std::max(q*Aref_*lRef_, tinyD);

    c_.Cd_p = fs*dot3(Fp_, e1_);
    c_.Cd_v = fs*dot3(Fv_, e1_);
    c_.Cl_p = fs*dot3(Fp_, e3_);
    c_.Cl_v = fs*dot3(Fv_, e3_);
    c_.Cd = c_.Cd_p + c_.Cd_v;
    c_.Cl = c_.Cl_p + c_.Cl_v;
    c_.Cm = ms*dot3(Mp_, e2_) + ms*dot3(Mv_, e2_);

    evaluated_ = true;
}


void Foam::coupledForces::evaluate(const volSymmTensorField& devReff)
{
    evaluate(mesh_.lookupObject<volScalarField>(pName_), devReff);
}


void Foam::coupledForces::write(const double iter)
{
    if (!evaluated_ || !UPstream::master())
    {
        return;
    }
    openFile();

    char buf[512];
    std::snprintf
    (
        buf, sizeof(buf),
        "%.10g %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
        iter, c_.Cd, c_.Cl, c_.Cm, c_.Cd_p, c_.Cd_v, c_.Cl_p, c_.Cl_v
    );
    std::ofstream os(file_.c_str(), std::ios::out | std::ios::app);
    os  << buf;
}


Foam::vector Foam::coupledForces::pressureForce() const
{
    return toVector(Fp_);
}


Foam::vector Foam::coupledForces::viscousForce() const
{
    return toVector(Fv_);
}


Foam::vector Foam::coupledForces::pressureMoment() const
{
    return toVector(Mp_);
}


Foam::vector Foam::coupledForces::viscousMoment() const
{
    return toVector(Mv_);
}


void Foam::coupledForces::writeSettings(dictionary& dict) const
{
    dictionary d;
    wordList names(patchIDs_.size());
    forAll(patchIDs_, i)
    {
        names[i] = mesh_.boundary()[patchIDs_[i]].name();
    }
    d.add("patches", names);
    d.add("rhoInf", rhoInf_);
    d.add("magUInf", magUInf_);
    d.add("lRef", lRef_);
    d.add("Aref", Aref_);
    d.add("pRef", pRef_);
    d.add("CofR", toVector(CofR_));
    d.add("e1", toVector(e1_));
    d.add("e2", toVector(e2_));
    d.add("e3", toVector(e3_));
    dict.add("coupledForces", d);
}


// ************************************************************************* //
