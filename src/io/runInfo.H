/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::runInfo

Description
    Once-per-run information (spec 12.2, 9.1):
    peak resident memory from /proc/self/status (VmHWM), flush-to-zero /
    denormals-are-zero control, FPE trap state, precision sizes.

SourceFiles
    runInfo.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_runInfo_H
#define coupledFoam_runInfo_H

#include "scalar.H"
#include "label.H"
#include "word.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace runInfo
{

//- Peak resident set size of this process in kB (VmHWM), 0 if unknown
label peakRSSkB();

//- CPU time (user + system) of this process in seconds
doubleScalar cpuSeconds();

//- Enable FTZ and DAZ (x86 SSE). Returns true if applied.
bool enableFTZ();

//- True if the FTZ bit is set in MXCSR
bool ftzActive();

//- True if FOAM_SIGFPE trapping is active
bool fpeActive();

//- Precision information as a dictionary
dictionary precisionDict();

} // End namespace runInfo
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
