#!/bin/bash
#------------------------------------------------------------------------------
# Test campaign D-063 (user, 2026-09-22) on one frozen commit.
#   setsid nohup /home/jonas/bin/cfenv sys bash bench/campaign.sh [phase...] &
# Phases (default: all, in this order; a failing phase is logged and the
# campaign continues):
#   freeze     tag + clean rebuild of the DP install and the SP install
#   archive    old run dirs/results aside (meshes and ref_* stay: cached refs)
#   battery    non-timing tests, concurrently in lanes on disjoint cores
#   timing     timing-relevant tests ALONE, sequentially: T0-T3 np1+np4, T4a
#   sp         SP (D-062/D-070): T0 Re100, T1 (gate override, informational),
#              T3-SST at 1 rank, T4a at 10 ranks - from an SP shell, alone
#   benchlight T1, T2, T3-SST, T3-GEKO: one lane per case, each on its own
#              physical core, nothing else running (disclosed in the paper)
#   bencht4a   T4a: A, B, C, H (1 repeat), alone
#   benchsp    F1/F2 (simpleFoam/coupledFoam SP): T1, T3-SST 3 repeats,
#              T4a 1 repeat, from an SP shell, alone
#   scaling    strong scaling on T4a (D-059), alone
#   t4b        only if the T4a test passed: T4b test + benchmark C, alone
#   report     make_report (guard mode, frozen commit) + both PDFs
# An interim report is built after timing, sp, benchlight and bencht4a so
# that a current PDF exists at any time (report_interim).
# No T5 (D-063). Log: run/campaign_<stamp>.log
#------------------------------------------------------------------------------
set -u
R=/home/jonas/coupledFoam
PY=/home/jonas/OF/venv/bin/python
STAMP=${CF_CAMPAIGN_STAMP:-$(date +%Y%m%d)}
LOG=$R/run/campaign_$STAMP.log
TAG=campaign-$STAMP
SPOF=/home/jonas/OpenFOAM-v2606-SP
PDFOUT=${CF_PDF_OUT:-/mnt/c/Users/DELLT5~1/AppData/Local/Temp/claude/C--Users-Dell-T5600-OneDrive-Dokumente-SIMS-Claude/4f47b3b5-3cfe-4a69-94f0-20066ca89e71/scratchpad}
# SP phases (sp, benchsp) are not in the default list: the SP simpleFoam
# references hit the float floor (2026-09-23 06:40, see D-070 follow-up)
DPCFG=A,B,C,D,F,G,H,H-tune,E-rcScalar,E-algPair,E-eta07,E-noSFD,E-nonOrth60,E-nonOrth65
phases=${*:-freeze archive battery timing benchlight bencht4a scaling t4b report}
say() { echo "[$(date '+%F %H:%M:%S')] $*" | tee -a $LOG; }
has() { case " $phases " in *" $1 "*) return 0;; esac; return 1; }
res() { grep -E "passed|failed|error|skipped|no tests ran" "$1" 2>/dev/null | tail -1; }
# a clean SP shell (cfenv refuses to stack environments)
sp() { env -i HOME=$HOME USER=$USER PATH=/usr/local/bin:/usr/bin:/bin \
           CFENV_CWD=$R "$@"; }
cd $R || exit 1
mkdir -p $R/run
SHA=$(git rev-parse HEAD)

report_interim() {   # report_interim <label>
    say "report ($1): guard mode, commit ${SHA:0:7}"
    nice -n 19 $PY bench/make_report.py --commit $SHA --allow-historical > $R/run/campaign_report_$1.log 2>&1
    local rr=$?
    (cd report/paper && nice -n 19 make > $R/run/campaign_make_$1.log 2>&1)
    local mr=$?
    cp report/paper/paper.pdf $PDFOUT/coupledFoam_paper_campaign.pdf 2>/dev/null
    cp report/paper/paper_tutorial.pdf $PDFOUT/coupledFoam_tutorial_campaign.pdf 2>/dev/null
    say "report ($1): make_report rc $rr ($(grep -ci pending $R/run/campaign_report_$1.log) pending notes), latex rc $mr, PDFs copied"
}

# --------------------------------------------------------------- freeze
if has freeze; then
    [ -z "$(git status --porcelain --untracked-files=no)" ] || { say "tree not clean - abort"; exit 1; }
    git tag -f $TAG > /dev/null
    say "freeze: $TAG = ${SHA:0:7}"
    (cd src && wclean > /dev/null 2>&1 && nice -n 19 wmake -j 12 libso > /tmp/cp_lib.log 2>&1) \
        || { say "DP lib build FAILED"; exit 1; }
    for a in applications/coupledFoam applications/test/Test-* applications/utilities/*; do
        [ -d "$a/Make" ] || continue
        (cd $a && wclean > /dev/null 2>&1; nice -n 19 wmake > /tmp/cp_app.log 2>&1) \
            || { say "DP build FAILED: $a"; exit 1; }
    done
    say "freeze: DP build ok"
    if has sp && [ "$(cat $SPOF/build_sp.done 2>/dev/null)" = "OK" ]; then
        sp CF_REPO=$R WM_NCOMPPROCS=12 /home/jonas/bin/cfenv sp bash $R/bench/build_sp_install.sh \
            > $R/run/campaign_sp_build.log 2>&1 \
            && say "freeze: SP build ok" || say "freeze: SP build FAILED (SP phases will fail)"
    else
        say "freeze: no SP OpenFOAM build"
    fi
fi

# -------------------------------------------------------------- archive
if has archive; then
    A=$R/run/_archive_$STAMP; mkdir -p $A
    for d in $R/run/*/; do
        n=$(basename $d)
        case $n in *_mesh*|ref_*|unit_cavity|_archive_*) continue;; esac
        mv "$d" "$A/"
    done
    RA=$R/results/_archive_$STAMP; mkdir -p $RA
    for sub in tests gates bench; do [ -d results/$sub ] && mv results/$sub $RA/$sub; done
    mkdir -p results/tests
    say "archive: run -> $A, results -> $RA"
fi

# -------------------------------------------------------------- battery
lane() {   # lane <name> <cpus> <pytest args...>
    local name=$1 cpus=$2; shift 2
    CF_MPI_BIND=none taskset -c $cpus nice -n 19 $PY -m pytest "$@" -q -p no:cacheprovider \
        > $R/run/campaign_$name.log 2>&1
    say "battery $name: $(res $R/run/campaign_$name.log)"
}
if has battery; then
    say "battery start (non-timing tests, lanes on disjoint physical cores)"
    lane unit     0,2,4,6     tests/test_unit.py tests/test_env.py --ranks 1,4 &
    lane gates    8           tests/test_gates.py tests/test_harness.py tests/test_keywords.py tests/test_coupledForces.py &
    lane restart  10,12,14,16 tests/test_restart.py --ranks 1,4 &
    lane fpe      18          tests/test_fpe.py --ranks 1 &
    lane diag     20,22,24,26 tests/test_diagnostics.py --ranks 1,4 &
    wait
    say "battery done"
fi

# --------------------------------------------------------------- timing
if has timing; then
    say "timing start (alone)"
    for t in test_T0_cavity.py test_T1_pitzDaily.py test_T2_backwardFacingStep.py test_T3_airFoil.py; do
        CF_MPI_BIND=core nice -n 19 $PY -m pytest tests/$t --ranks 1,4 -q -p no:cacheprovider \
            > $R/run/campaign_timing_${t%.py}.log 2>&1
        say "timing $t: $(res $R/run/campaign_timing_${t%.py}.log)"
    done
    CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 nice -n 19 $PY -m pytest \
        "tests/test_T4_motorBike.py::test_T4[a]" --heavy -q -p no:cacheprovider \
        > $R/run/campaign_timing_T4a.log 2>&1
    say "timing T4a: $(res $R/run/campaign_timing_T4a.log)"
    report_interim timing
fi

# ------------------------------------------------------------------- sp
spt() {   # spt <log> <pytest args...>   - pytest in a clean SP shell
    local lg=$1; shift
    sp CF_PRECISION=sp CF_MPI_BIND=core "$@" > $R/run/campaign_sp_$lg.log 2>&1
    say "sp $lg: $(res $R/run/campaign_sp_$lg.log)"
}
if has sp; then
    say "sp start (alone; D-070 settings)"
    spt T0 /home/jonas/bin/cfenv sp $PY -m pytest tests/test_T0_cavity.py --ranks 1 -k Re100 -q -p no:cacheprovider
    spt T1 env CF_SP_GATE_OVERRIDE="Boundary openness" /home/jonas/bin/cfenv sp $PY -m pytest \
        tests/test_T1_pitzDaily.py --ranks 1 -q -p no:cacheprovider
    spt T3 /home/jonas/bin/cfenv sp $PY -m pytest "tests/test_T3_airFoil.py::test_T3[np1-kOmegaSST]" \
        --ranks 1 -q -p no:cacheprovider
    spt T4a env CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 /home/jonas/bin/cfenv sp $PY -m pytest \
        "tests/test_T4_motorBike.py::test_T4[a]" --heavy -q -p no:cacheprovider
    report_interim sp
fi

# ----------------------------------------------------------- benchlight
if has benchlight; then
    say "benchlight start: one lane per case, each on its own physical core, nothing else running"
    i=0
    for c in T1 T2 T3-SST T3-GEKO; do
        cpu=$((2 * (12 + i)))      # physical cores 12..15, one logical CPU each
        ( CF_MPI_BIND=none taskset -c $cpu nice -n 19 $PY bench/run_bench.py --cases $c --configs $DPCFG --repeats 3 --allow-busy \
              > $R/run/campaign_bench_$c.log 2>&1
          say "benchlight $c rc $?" ) &
        i=$((i + 1))
    done
    wait
    say "benchlight done"
    report_interim benchlight
fi

# ------------------------------------------------------------- bencht4a
if has bencht4a; then
    say "bencht4a start (alone)"
    CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 nice -n 19 $PY bench/run_bench.py --cases T4a --configs $DPCFG --repeats 1 \
        > $R/run/campaign_bench_T4a.log 2>&1
    say "bencht4a rc $?"
    report_interim bencht4a
fi

# -------------------------------------------------------------- benchsp
if has benchsp; then
    say "benchsp start (SP shell, alone)"
    sp CF_PRECISION=sp CF_MPI_BIND=none /home/jonas/bin/cfenv sp $PY bench/run_bench.py \
        --cases T1,T3-SST --configs F1,F2 --repeats 3 > $R/run/campaign_benchsp_light.log 2>&1
    say "benchsp light rc $?"
    sp CF_PRECISION=sp CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 /home/jonas/bin/cfenv sp $PY bench/run_bench.py \
        --cases T4a --configs F1,F2 --repeats 1 > $R/run/campaign_benchsp_T4a.log 2>&1
    say "benchsp T4a rc $?"
fi

# -------------------------------------------------------------- scaling
if has scaling; then
    say "scaling start (alone, T4a, D-059)"
    CF_SCALING_RANKS=1,2,4,8,12,16 CF_SCALING_MESH=a CF_SCALING_ITERS=150 \
    CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_BIND=core nice -n 19 \
        $PY -m pytest tests/test_scaling.py --heavy -q -p no:cacheprovider > $R/run/campaign_scaling.log 2>&1
    say "scaling: $(res $R/run/campaign_scaling.log)"
fi

# ------------------------------------------------------------------ t4b
if has t4b; then
    if grep -qE "1 passed" $R/run/campaign_timing_T4a.log 2>/dev/null; then
        say "t4b start (T4a passed; alone)"
        CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 nice -n 19 $PY -m pytest \
            "tests/test_T4_motorBike.py::test_T4[b]" --heavy -q -p no:cacheprovider \
            > $R/run/campaign_timing_T4b.log 2>&1
        say "t4b test: $(res $R/run/campaign_timing_T4b.log)"
        CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 nice -n 19 $PY bench/run_bench.py --cases T4b --configs C --repeats 1 \
            > $R/run/campaign_bench_T4b.log 2>&1
        say "t4b bench rc $?"
    else
        say "t4b skipped: T4a test did not pass ($(res $R/run/campaign_timing_T4a.log))"
    fi
fi

# --------------------------------------------------------------- report
if has report; then
    report_interim final
fi
say "campaign finished: $phases"
