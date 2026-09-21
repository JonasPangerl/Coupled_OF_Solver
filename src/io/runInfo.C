/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "runInfo.H"
#include "blockScalar.H"
#include "sigFpe.H"
#include <fstream>
#include <string>
#include <cstdlib>
#include <sys/resource.h>

#if defined(__SSE__) || defined(__x86_64__)
    #include <xmmintrin.h>
    #include <pmmintrin.h>
    #define COUPLEDFOAM_HAVE_SSE 1
#endif

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

Foam::label Foam::runInfo::peakRSSkB()
{
    std::ifstream is("/proc/self/status");
    std::string line;
    while (std::getline(is, line))
    {
        if (line.compare(0, 6, "VmHWM:") == 0)
        {
            return label(std::strtol(line.c_str() + 6, nullptr, 10));
        }
    }
    return 0;
}


double Foam::runInfo::cpuSeconds()
{
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0)
    {
        return 0;
    }
    return
        double(ru.ru_utime.tv_sec) + 1e-6*double(ru.ru_utime.tv_usec)
      + double(ru.ru_stime.tv_sec) + 1e-6*double(ru.ru_stime.tv_usec);
}


bool Foam::runInfo::enableFTZ()
{
#ifdef COUPLEDFOAM_HAVE_SSE
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    return true;
#else
    return false;
#endif
}


bool Foam::runInfo::ftzActive()
{
#ifdef COUPLEDFOAM_HAVE_SSE
    return _MM_GET_FLUSH_ZERO_MODE() == _MM_FLUSH_ZERO_ON;
#else
    return false;
#endif
}


bool Foam::runInfo::fpeActive()
{
    return sigFpe::active();
}


Foam::dictionary Foam::runInfo::precisionDict()
{
    dictionary d;
    d.add("sizeofScalar", label(sizeof(scalar)));
    d.add("sizeofSolveScalar", label(sizeof(solveScalar)));
    d.add("sizeofBlockScalar", label(sizeof(blockScalar)));
    d.add("sizeofReduceScalar", label(sizeof(reduceScalar)));
    const char* opts = std::getenv("WM_OPTIONS");
    d.add("WM_OPTIONS", word(opts ? opts : "unset"));
    return d;
}


// ************************************************************************* //
