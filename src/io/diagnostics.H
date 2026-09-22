/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::diagJson
    Foam::diagnostics

Description
    Staged deep-diagnostics logging (OPUS_TASKS TASK 5, DECISIONS D-045).

    Dictionary coupled.diagnostics:
        level     0 | 1 | 2 | 3   (default 0: off)
        echo      yes | no        (default no: print the level-1 record
                                   pretty-printed into the solver log)
        maxBytes  bytes per rank  (default 2 GiB; when exceeded one final
                                   truncation record is written and logging
                                   stops, the run continues)
        upLeg     first | all     (default first: level-2 up-leg GAMG norms
                                   on the first preconditioner application
                                   of each solve only; level 3: all)

    Output (level >= 1): <case>/diagnostics/diag.rank<N>.jsonl, one JSON
    object per line (JSON Lines), numbers with max_digits10, flushed once per
    outer iteration. Line 1 of every file is a header record
    ("type": "header"), every outer iteration writes one "type": "iter"
    record, a size overflow one "type": "truncated" record.

    Levels
        1  one record per outer iteration (residuals, controls, turbulence,
           timings, GAMG hierarchy, memory, phase)
        2  + per linear solve: Krylov residual history, true residual at
           every restart, per preconditioner application the per-level
           residual norms of the down leg (before/after pre-smoothing) and,
           for the first application of each solve (upLeg first, default)
           or all (upLeg all), of the up leg (before/after post-smoothing),
           coarsest-level its/residual, K-cycle GCR values, rhoOpt of the
           first application
        3  + per smoother sweep residual norms, per-level operator stats on
           set-up, line-search trial violations, per-patch inflow/outflow
           flips, Anderson internals

    HARD RULE (TASK 5.1): level 0 builds no strings, gathers no fields and
    adds no collectives. Every hook is guarded at the call site with
    active(n); the modules receive a pointer that is null (or inactive) at
    level 0.

    diagJson is a minimal streaming JSON builder (no DOM): containers are
    opened and closed explicitly, commas are inserted automatically.

SourceFiles
    diagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_diagnostics_H
#define coupledFoam_diagnostics_H

#include "label.H"
#include "scalar.H"
#include "dictionary.H"
#include "fileName.H"
#include "DynamicList.H"
#include <string>
#include <fstream>
#include <chrono>
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                          Class diagJson Declaration
\*---------------------------------------------------------------------------*/

class diagJson
{
    // Private Data

        std::string s_;

        //- One entry per open container: true while it is still empty
        DynamicList<bool> empty_;


    // Private Member Functions

        //- Comma before a new element if the container is not empty
        void sep();

        void key(const char* k);


public:

    // Static Functions

        //- Append a double with max_digits10 (17), null if not finite
        static void appendNumber(std::string& s, const doubleScalar v);

        //- Append a quoted, escaped string
        static void appendString(std::string& s, const char* v);

        //- Re-indent a compact JSON object for the log (objects broken
        //  over lines, arrays kept on one line)
        static std::string pretty(const std::string& json);


    // Constructors

        diagJson();


    // Member Functions

        void clear();

        //- Containers (keyed: inside an object; unkeyed: inside an array
        //  or top level)
        void beginObject();
        void beginObject(const char* k);
        void endObject();
        void beginArray(const char* k);
        void beginArray();
        void endArray();

        //- Keyed values
        void add(const char* k, const doubleScalar v);
        void add(const char* k, const label v);
        void add(const char* k, const bool v);
        void add(const char* k, const char* v);
        void add(const char* k, const std::string& v);
        void addNull(const char* k);

        //- Keyed, already serialised JSON value
        void addRaw(const char* k, const std::string& json);

        //- Array elements
        void value(const doubleScalar v);
        void value(const label v);
        void valueRaw(const std::string& json);

        //- Keyed arrays
        void addList(const char* k, const UList<doubleScalar>& v);
        void addList(const char* k, const UList<label>& v);

        const std::string& str() const noexcept { return s_; }
        bool closed() const noexcept { return empty_.empty(); }
};


/*---------------------------------------------------------------------------*\
                         Class diagnostics Declaration
\*---------------------------------------------------------------------------*/

class diagnostics
{
    // Private Data

        label level_;
        bool echo_;
        std::int64_t maxBytes_;

        //- Up-leg level norms on every application (else the first of
        //  each solve only; level 3: always every application)
        bool upLegAll_;

        //- Output file of this rank
        fileName file_;
        std::ofstream os_;
        std::int64_t bytes_;
        bool opened_;
        bool truncated_;

        // Level 1: per-iteration accumulators

            doubleScalar tPrecSetup_;
            doubleScalar tPrecApply_;
            label nPrecApply_;
            label nPrecSetup_;
            doubleScalar tDiag_;

        // Level 2: linear-solve log of the current outer iteration

            //- Solve objects of this iteration, comma separated
            std::string solves_;
            label nSolves_;

            //- Current solve
            DynamicList<doubleScalar> krylov_;
            DynamicList<doubleScalar> restartRes_;
            doubleScalar rhoOpt_;
            std::string apps_;
            label nApps_;

            //- Current preconditioner application
            std::string appLevels_;
            std::string appCoarse_;
            std::string appK_;

        // Level 3

            //- Per-level operator stats of set-ups in this iteration
            std::string setup_;

        //- Extra (level >= 2) keyed fields of the current record, added by
        //  the solver application (object body without braces)
        diagJson ext_;


    // Private Member Functions

        //- Write one line (bytes accounting, truncation)
        void writeLine(const std::string& line);


public:

    // Static Functions

        //- Monotonic wall clock in seconds (cheap, vDSO)
        static doubleScalar clock()
        {
            return std::chrono::duration<doubleScalar>
            (
                std::chrono::steady_clock::now().time_since_epoch()
            ).count();
        }

        //- Current resident set size of this process [kB] (/proc/self/statm)
        static label currentRSSkB();

        //- Peak resident set size of this process [kB] (getrusage)
        static label peakRSSkB();


    // Constructors

        //- Construct from the "coupled" dictionary (subdictionary
        //  diagnostics); does not open the file
        explicit diagnostics(const dictionary& coupledDict);

        //- No copy
        diagnostics(const diagnostics&) = delete;
        void operator=(const diagnostics&) = delete;


    // Member Functions

        label level() const noexcept { return level_; }
        bool active(const label l) const noexcept { return level_ >= l; }
        bool echo() const noexcept { return echo_; }

        //- Log the up-leg (post-smoothing) level norms in the current
        //  preconditioner application
        bool upLeg() const noexcept
        {
            return level_ >= 3 || upLegAll_ || nApps_ == 0;
        }
        bool truncated() const noexcept { return truncated_; }
        const fileName& file() const noexcept { return file_; }

        //- Level >= 1 only: create <case>/diagnostics and open the rank
        //  file (append: continue a restarted run), write the header
        //  record (keyed fields given as an already built object)
        void open
        (
            const Time& runTime,
            const bool append,
            const std::string& headerJson
        );

        //- Reset the per-iteration accumulators and buffers
        void beginIteration();


        // Level 1 accumulators (preconditioner timing, diagnostics cost)

            void addPrecSetup(const doubleScalar t)
            {
                tPrecSetup_ += t;
                ++nPrecSetup_;
            }
            void addPrecApply(const doubleScalar t)
            {
                tPrecApply_ += t;
                ++nPrecApply_;
            }
            void addDiagTime(const doubleScalar t) { tDiag_ += t; }

            doubleScalar tPrecSetup() const noexcept { return tPrecSetup_; }
            doubleScalar tPrecApply() const noexcept { return tPrecApply_; }
            label nPrecApply() const noexcept { return nPrecApply_; }
            label nPrecSetup() const noexcept { return nPrecSetup_; }
            doubleScalar tDiag() const noexcept { return tDiag_; }


        // Level 2 hooks (linear solver, blockGAMG)

            void beginSolve();
            void krylovResidual(const doubleScalar r);
            void restartResidual(const doubleScalar r);
            void rhoOpt(const doubleScalar r);

            //- Close the current solve with its performance values
            void endSolve
            (
                const char* solverName,
                const doubleScalar initialResidual,
                const doubleScalar finalResidual,
                const label nIterations,
                const label nRestarts,
                const bool converged,
                const doubleScalar rho
            );

            void beginApplication();

            //- One visit of level l in a cycle: residual norm before and
            //  after pre-smoothing (down leg), before and after
            //  post-smoothing (up leg, if up: see upLeg()); per-sweep norms
            //  (level 3, may be empty)
            void levelVisit
            (
                const label l,
                const doubleScalar preBefore,
                const doubleScalar preAfter,
                const bool up,
                const doubleScalar postBefore,
                const doubleScalar postAfter,
                const UList<doubleScalar>& preSweeps,
                const UList<doubleScalar>& postSweeps
            );

            //- Coarsest-level solve (dense: direct LU, its 0, no residual)
            void coarseSolve
            (
                const label its,
                const doubleScalar finalResidual,
                const bool dense
            );

            //- K-cycle GCR step at level l
            void kStep
            (
                const label l,
                const doubleScalar r0,
                const doubleScalar r1,
                const doubleScalar threshold,
                const bool secondStep,
                const doubleScalar a1,
                const doubleScalar a2
            );

            void endApplication();


        // Level 3 hooks

            //- Operator statistics of level l on a set-up
            void operatorStats
            (
                const label l,
                const label rows,
                const doubleScalar dominanceMin,
                const doubleScalar dominanceMedian,
                const label nRowsNoOffDiag
            );


        // Record

            //- Keyed fields added to the current record at level >= 2
            //  (object body; the solver application adds level-3 items)
            diagJson& ext() noexcept { return ext_; }

            //- Finish and write the record: level1 is a closed JSON object
            //  holding the level-1 fields; the level-2/3 buffers are
            //  appended. echo prints level1 pretty-printed (master).
            void writeRecord(const std::string& level1);


        //- Serialise the linear-solve buffers of this iteration (level 2)
        std::string linearJson() const;
};


/*---------------------------------------------------------------------------*\
                       Class diagPhase Declaration
\*---------------------------------------------------------------------------*/

//- Phase classification of TASK 5.2, rules in this order:
//    1. startup     beta < 1
//    2. stalled     no new minimum of R in the last stallWindow iterations
//                   (min R of the window >= min R of everything before it)
//    3. asymptotic  R <= resFactor*residualTol (residualTol > 0), or the
//                   Cd/Cl window criterion within forceFactor of its tolerance
//    4. ramp        otherwise
class diagPhase
{
    label window_;
    doubleScalar resFactor_;
    doubleScalar forceFactor_;

    //- R of the last window_ iterations (ring buffer)
    DynamicList<doubleScalar> ring_;
    label head_;

    //- min R of all iterations older than the window
    doubleScalar bestBefore_;

public:

    diagPhase();

    //- Record R of this iteration and classify. forceRatio: spread of
    //  the Cd/Cl window over its tolerance (< 0: not available).
    const char* classify
    (
        const doubleScalar beta,
        const doubleScalar R,
        const doubleScalar residualTol,
        const doubleScalar forceRatio
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
