/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "diagnostics.H"
#include "coupledDefaults.H"
#include "Time.H"
#include "OSspecific.H"
#include "UPstream.H"
#include "IOstreams.H"
#include <cmath>
#include <cstdio>
#include <vector>
#include <sys/resource.h>
#include <unistd.h>

// * * * * * * * * * * * * * * * * diagJson  * * * * * * * * * * * * * * * * //

void Foam::diagJson::appendNumber(std::string& s, const doubleScalar v)
{
    if (!std::isfinite(v))
    {
        s += "null";
        return;
    }
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    s += buf;
}


void Foam::diagJson::appendString(std::string& s, const char* v)
{
    s += '"';
    for (const char* p = v; *p; ++p)
    {
        switch (*p)
        {
            case '"':  s += "\\\""; break;
            case '\\': s += "\\\\"; break;
            case '\n': s += "\\n"; break;
            case '\t': s += "\\t"; break;
            default:   s += *p; break;
        }
    }
    s += '"';
}


std::string Foam::diagJson::pretty(const std::string& json)
{
    // Objects outside any array are broken over lines, everything inside
    // an array stays inline
    const std::string ind("  ");
    std::string out;
    std::vector<char> stack;
    label nArrays = 0;
    label depth = 0;
    bool inStr = false;

    auto newline = [&](const label d)
    {
        out += '\n';
        for (label i = 0; i < d; ++i)
        {
            out += ind;
        }
    };

    for (std::size_t i = 0; i < json.size(); ++i)
    {
        const char c = json[i];
        if (inStr)
        {
            out += c;
            if (c == '\\' && i + 1 < json.size())
            {
                out += json[++i];
            }
            else if (c == '"')
            {
                inStr = false;
            }
            continue;
        }
        switch (c)
        {
            case '"':
                inStr = true;
                out += c;
                break;
            case '{':
                out += c;
                stack.push_back(c);
                if (nArrays == 0)
                {
                    newline(++depth);
                }
                break;
            case '[':
                out += c;
                stack.push_back(c);
                ++nArrays;
                break;
            case '}':
                if (!stack.empty())
                {
                    stack.pop_back();
                }
                if (nArrays == 0)
                {
                    newline(--depth);
                }
                out += c;
                break;
            case ']':
                if (!stack.empty())
                {
                    stack.pop_back();
                }
                --nArrays;
                out += c;
                break;
            case ',':
                out += c;
                if (nArrays == 0)
                {
                    newline(depth);
                }
                else
                {
                    out += ' ';
                }
                break;
            case ':':
                out += ": ";
                break;
            default:
                out += c;
                break;
        }
    }
    return out;
}


Foam::diagJson::diagJson()
:
    s_(),
    empty_()
{}


void Foam::diagJson::clear()
{
    s_.clear();
    empty_.clear();
}


void Foam::diagJson::sep()
{
    if (empty_.size())
    {
        if (!empty_.last())
        {
            s_ += ',';
        }
        empty_.last() = false;
    }
}


void Foam::diagJson::key(const char* k)
{
    sep();
    appendString(s_, k);
    s_ += ':';
}


void Foam::diagJson::beginObject()
{
    sep();
    s_ += '{';
    empty_.append(true);
}


void Foam::diagJson::beginObject(const char* k)
{
    key(k);
    s_ += '{';
    empty_.append(true);
}


void Foam::diagJson::endObject()
{
    s_ += '}';
    empty_.remove();
}


void Foam::diagJson::beginArray(const char* k)
{
    key(k);
    s_ += '[';
    empty_.append(true);
}


void Foam::diagJson::beginArray()
{
    sep();
    s_ += '[';
    empty_.append(true);
}


void Foam::diagJson::endArray()
{
    s_ += ']';
    empty_.remove();
}


void Foam::diagJson::add(const char* k, const doubleScalar v)
{
    key(k);
    appendNumber(s_, v);
}


void Foam::diagJson::add(const char* k, const label v)
{
    key(k);
    s_ += std::to_string(v);
}


void Foam::diagJson::add(const char* k, const bool v)
{
    key(k);
    s_ += (v ? "true" : "false");
}


void Foam::diagJson::add(const char* k, const char* v)
{
    key(k);
    appendString(s_, v);
}


void Foam::diagJson::add(const char* k, const std::string& v)
{
    add(k, v.c_str());
}


void Foam::diagJson::addNull(const char* k)
{
    key(k);
    s_ += "null";
}


void Foam::diagJson::addRaw(const char* k, const std::string& json)
{
    key(k);
    s_ += json;
}


void Foam::diagJson::value(const doubleScalar v)
{
    sep();
    appendNumber(s_, v);
}


void Foam::diagJson::value(const label v)
{
    sep();
    s_ += std::to_string(v);
}


void Foam::diagJson::valueRaw(const std::string& json)
{
    sep();
    s_ += json;
}


void Foam::diagJson::addList(const char* k, const UList<doubleScalar>& v)
{
    beginArray(k);
    for (const doubleScalar x : v)
    {
        value(x);
    }
    endArray();
}


void Foam::diagJson::addList(const char* k, const UList<label>& v)
{
    beginArray(k);
    for (const label x : v)
    {
        value(x);
    }
    endArray();
}


// * * * * * * * * * * * * * * * * diagnostics * * * * * * * * * * * * * * * //

Foam::label Foam::diagnostics::currentRSSkB()
{
    long pages = 0;
    long resident = 0;
    std::FILE* f = std::fopen("/proc/self/statm", "r");
    if (!f)
    {
        return 0;
    }
    const int n = std::fscanf(f, "%ld %ld", &pages, &resident);
    std::fclose(f);
    if (n != 2)
    {
        return 0;
    }
    const long pageKB = sysconf(_SC_PAGESIZE)/1024;
    return label(resident*pageKB);
}


Foam::label Foam::diagnostics::peakRSSkB()
{
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0)
    {
        return 0;
    }
    // Linux: ru_maxrss in kB
    return label(ru.ru_maxrss);
}


Foam::diagnostics::diagnostics(const dictionary& coupledDict)
:
    level_(coupledDefaults::diagLevel),
    echo_(coupledDefaults::diagEcho),
    maxBytes_(coupledDefaults::diagMaxBytes),
    file_(),
    os_(),
    bytes_(0),
    opened_(false),
    truncated_(false),
    tPrecSetup_(0),
    tPrecApply_(0),
    nPrecApply_(0),
    nPrecSetup_(0),
    tDiag_(0),
    solves_(),
    nSolves_(0),
    krylov_(),
    restartRes_(),
    rhoOpt_(-1),
    apps_(),
    nApps_(0),
    appLevels_(),
    appCoarse_(),
    appK_(),
    setup_(),
    ext_()
{
    const dictionary& d = coupledDict.subOrEmptyDict("diagnostics");
    level_ = d.getOrDefault<label>("level", coupledDefaults::diagLevel);
    echo_ = d.getOrDefault<bool>("echo", coupledDefaults::diagEcho);
    maxBytes_ =
        d.getOrDefault<std::int64_t>("maxBytes", coupledDefaults::diagMaxBytes);

    if (level_ < 0 || level_ > coupledDefaults::diagMaxLevel)
    {
        FatalIOErrorInFunction(d)
            << "coupled.diagnostics.level must be 0.."
            << coupledDefaults::diagMaxLevel << ", got " << level_
            << exit(FatalIOError);
    }
    if (maxBytes_ <= 0)
    {
        FatalIOErrorInFunction(d)
            << "coupled.diagnostics.maxBytes must be > 0"
            << exit(FatalIOError);
    }
}


void Foam::diagnostics::open
(
    const Time& runTime,
    const bool append,
    const std::string& headerJson
)
{
    if (!active(1) || opened_)
    {
        return;
    }

    const fileName dir(runTime.globalPath()/"diagnostics");
    Foam::mkDir(dir);
    file_ = dir/("diag.rank" + Foam::name(UPstream::myProcNo()) + ".jsonl");

    os_.open
    (
        file_,
        append ? (std::ios::out | std::ios::app) : std::ios::out
    );
    if (!os_.good())
    {
        FatalErrorInFunction
            << "cannot open diagnostics file " << file_ << exit(FatalError);
    }
    opened_ = true;
    bytes_ = (append ? std::int64_t(os_.tellp()) : 0);
    if (bytes_ < 0)
    {
        bytes_ = 0;
    }

    writeLine(headerJson);
}


void Foam::diagnostics::writeLine(const std::string& line)
{
    if (!opened_ || truncated_)
    {
        return;
    }

    const std::int64_t n = std::int64_t(line.size()) + 1;
    if (bytes_ + n > maxBytes_)
    {
        // One final truncation record (not counted against the limit),
        // then stop logging; the run continues
        std::string t("{\"type\":\"truncated\",\"bytesWritten\":");
        t += std::to_string(bytes_);
        t += ",\"maxBytes\":";
        t += std::to_string(maxBytes_);
        t += ",\"rejectedRecordBytes\":";
        t += std::to_string(n);
        t += '}';
        os_ << t << '\n';
        os_.flush();
        truncated_ = true;
        return;
    }

    os_ << line << '\n';
    os_.flush();
    bytes_ += n;
}


void Foam::diagnostics::beginIteration()
{
    tPrecSetup_ = 0;
    tPrecApply_ = 0;
    nPrecApply_ = 0;
    nPrecSetup_ = 0;
    tDiag_ = 0;
    solves_.clear();
    nSolves_ = 0;
    setup_.clear();
    ext_.clear();
    ext_.beginObject();
}


// * * * * * * * * * * * * * * * * Level 2 * * * * * * * * * * * * * * * * * //

void Foam::diagnostics::beginSolve()
{
    krylov_.clear();
    restartRes_.clear();
    rhoOpt_ = -1;
    apps_.clear();
    nApps_ = 0;
}


void Foam::diagnostics::krylovResidual(const doubleScalar r)
{
    krylov_.append(r);
}


void Foam::diagnostics::restartResidual(const doubleScalar r)
{
    restartRes_.append(r);
}


void Foam::diagnostics::rhoOpt(const doubleScalar r)
{
    rhoOpt_ = r;
}


void Foam::diagnostics::endSolve
(
    const char* solverName,
    const doubleScalar initialResidual,
    const doubleScalar finalResidual,
    const label nIterations,
    const label nRestarts,
    const bool converged,
    const doubleScalar rho
)
{
    const doubleScalar t0 = clock();

    diagJson j;
    j.beginObject();
    j.add("solver", solverName);
    j.add("initial", initialResidual);
    j.add("final", finalResidual);
    j.add("its", nIterations);
    j.add("restarts", nRestarts);
    j.add("converged", converged);
    j.add("rho", rho);
    if (rhoOpt_ >= 0)
    {
        j.add("rhoOpt", rhoOpt_);
    }
    else
    {
        j.addNull("rhoOpt");
    }
    j.addList("krylov", krylov_);
    j.addList("trueResidualAtRestart", restartRes_);
    j.addRaw("precon", "[" + apps_ + "]");
    j.endObject();

    if (nSolves_++)
    {
        solves_ += ',';
    }
    solves_ += j.str();

    tDiag_ += clock() - t0;
}


void Foam::diagnostics::beginApplication()
{
    appLevels_.clear();
    appCoarse_.clear();
    appK_.clear();
}


void Foam::diagnostics::levelVisit
(
    const label l,
    const doubleScalar preBefore,
    const doubleScalar preAfter,
    const doubleScalar postBefore,
    const doubleScalar postAfter,
    const UList<doubleScalar>& preSweeps,
    const UList<doubleScalar>& postSweeps
)
{
    std::string& s = appLevels_;
    if (!s.empty())
    {
        s += ',';
    }
    s += "{\"l\":";
    s += std::to_string(l);
    s += ",\"pre\":[";
    diagJson::appendNumber(s, preBefore);
    s += ',';
    diagJson::appendNumber(s, preAfter);
    s += "],\"post\":[";
    diagJson::appendNumber(s, postBefore);
    s += ',';
    diagJson::appendNumber(s, postAfter);
    s += ']';

    auto list = [&s](const char* k, const UList<doubleScalar>& v)
    {
        if (v.empty())
        {
            return;
        }
        s += ",\"";
        s += k;
        s += "\":[";
        forAll(v, i)
        {
            if (i)
            {
                s += ',';
            }
            diagJson::appendNumber(s, v[i]);
        }
        s += ']';
    };
    list("preSweeps", preSweeps);
    list("postSweeps", postSweeps);
    s += '}';
}


void Foam::diagnostics::coarseSolve
(
    const label its,
    const doubleScalar finalResidual,
    const bool dense
)
{
    std::string& s = appCoarse_;
    if (!s.empty())
    {
        s += ',';
    }
    if (dense)
    {
        s += "{\"dense\":true,\"its\":0}";
        return;
    }
    s += "{\"its\":";
    s += std::to_string(its);
    s += ",\"res\":";
    diagJson::appendNumber(s, finalResidual);
    s += '}';
}


void Foam::diagnostics::kStep
(
    const label l,
    const doubleScalar r0,
    const doubleScalar r1,
    const doubleScalar threshold,
    const bool secondStep,
    const doubleScalar a1,
    const doubleScalar a2
)
{
    std::string& s = appK_;
    if (!s.empty())
    {
        s += ',';
    }
    s += "{\"l\":";
    s += std::to_string(l);
    s += ",\"r0\":";
    diagJson::appendNumber(s, r0);
    s += ",\"r1\":";
    diagJson::appendNumber(s, r1);
    s += ",\"threshold\":";
    diagJson::appendNumber(s, threshold);
    s += ",\"second\":";
    s += (secondStep ? "true" : "false");
    s += ",\"a1\":";
    diagJson::appendNumber(s, a1);
    s += ",\"a2\":";
    if (secondStep)
    {
        diagJson::appendNumber(s, a2);
    }
    else
    {
        s += "null";
    }
    s += '}';
}


void Foam::diagnostics::endApplication()
{
    if (nApps_++)
    {
        apps_ += ',';
    }
    apps_ += "{\"levels\":[";
    apps_ += appLevels_;
    apps_ += "],\"coarse\":[";
    apps_ += appCoarse_;
    apps_ += "],\"k\":[";
    apps_ += appK_;
    apps_ += "]}";
}


// * * * * * * * * * * * * * * * * Level 3 * * * * * * * * * * * * * * * * * //

void Foam::diagnostics::operatorStats
(
    const label l,
    const label rows,
    const doubleScalar dominanceMin,
    const doubleScalar dominanceMedian,
    const label nRowsNoOffDiag
)
{
    std::string& s = setup_;
    if (!s.empty())
    {
        s += ',';
    }
    s += "{\"l\":";
    s += std::to_string(l);
    s += ",\"rows\":";
    s += std::to_string(rows);
    s += ",\"dominanceMin\":";
    diagJson::appendNumber(s, dominanceMin);
    s += ",\"dominanceMedian\":";
    diagJson::appendNumber(s, dominanceMedian);
    s += ",\"rowsWithoutOffDiagonal\":";
    s += std::to_string(nRowsNoOffDiag);
    s += '}';
}


// * * * * * * * * * * * * * * * * * Record  * * * * * * * * * * * * * * * * //

std::string Foam::diagnostics::linearJson() const
{
    return "[" + solves_ + "]";
}


void Foam::diagnostics::writeRecord(const std::string& level1)
{
    if (!active(1))
    {
        return;
    }

    if (echo_ && UPstream::master())
    {
        Info<< "CFdiag| " << diagJson::pretty(level1).c_str() << endl;
    }

    if (!opened_ || truncated_)
    {
        return;
    }

    // level1 is a closed object: drop the '}' and append the other levels
    std::string rec(level1, 0, level1.size() - 1);

    if (active(2))
    {
        rec += ",\"linear\":";
        rec += linearJson();
    }
    if (active(3))
    {
        rec += ",\"gamgSetup\":[";
        rec += setup_;
        rec += ']';
    }
    ext_.endObject();
    const std::string& e = ext_.str();
    if (e.size() > 2)
    {
        rec += ',';
        rec.append(e, 1, e.size() - 2);
    }
    rec += '}';

    writeLine(rec);
}


// * * * * * * * * * * * * * * * * diagPhase * * * * * * * * * * * * * * * * //

Foam::diagPhase::diagPhase()
:
    window_(coupledDefaults::diagStallWindow),
    resFactor_(coupledDefaults::diagAsymptoticResidualFactor),
    forceFactor_(coupledDefaults::diagAsymptoticForceFactor),
    ring_(),
    head_(0),
    bestBefore_(GREAT)
{}


const char* Foam::diagPhase::classify
(
    const doubleScalar beta,
    const doubleScalar R,
    const doubleScalar residualTol,
    const doubleScalar forceRatio
)
{
    // History first: every iteration counts, whatever its phase
    if (ring_.size() < window_)
    {
        ring_.append(R);
    }
    else
    {
        bestBefore_ = min(bestBefore_, ring_[head_]);
        ring_[head_] = R;
        head_ = (head_ + 1) % window_;
    }

    // 1. start-up (upwind / deferred-correction ramp)
    if (beta < 1)
    {
        return "startup";
    }

    // 2. stalled: no new minimum within the last window iterations
    if (ring_.size() == window_ && bestBefore_ < GREAT)
    {
        doubleScalar wmin = GREAT;
        for (const doubleScalar r : ring_)
        {
            wmin = min(wmin, r);
        }
        if (wmin >= bestBefore_)
        {
            return "stalled";
        }
    }

    // 3. asymptotic
    if
    (
        (residualTol > 0 && R <= resFactor_*residualTol)
     || (forceRatio >= 0 && forceRatio <= forceFactor_)
    )
    {
        return "asymptotic";
    }

    // 4. ramp
    return "ramp";
}


// ************************************************************************* //
