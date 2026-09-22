#!/bin/bash
#------------------------------------------------------------------------------
# Final clean re-run (OPUS_TASKS TASK 6). Run from WSL:
#   /home/jonas/bin/cfenv sys bash final_rerun.sh <phase...>
# phases: freeze archive heavy light scaling bench report   (default: all)
#
# Lanes (16 physical cores, SMT siblings adjacent: logical 2k,2k+1 = core k):
#   heavy lane: T4a, T4b, T5 sequentially, 10 ranks bound to physical cores 0-9
#               (mpirun --bind-to core --map-by core, timing-clean)
#   light lane: T0-T3, restart, fpe, diagnostics, unit battery, pinned with
#               taskset to one logical CPU on each of physical cores 10-15
#               (logical 20 22 24 26 28 30), i.e. no SMT sharing with the
#               heavy lane or among themselves
# Every result is produced by the frozen commit; the report runs in guard mode.
#------------------------------------------------------------------------------
set -u
R=/home/jonas/coupledFoam
PY=/home/jonas/OF/venv/bin/python
STAMP=$(date +%Y%m%d)
LOG=$R/run/final_rerun_$STAMP.log
LIGHT_CPUS=20,22,24,26,28,30
phases=${*:-freeze archive heavy light scaling bench report}
say() { echo "[$(date '+%H:%M:%S')] $*" | tee -a $LOG; }
cd $R || exit 1

has() { case " $phases " in *" $1 "*) return 0;; esac; return 1; }

if has freeze; then
    [ -z "$(git status --porcelain --untracked-files=no)" ] || { say "tree not clean - abort"; exit 1; }
    TAG=final-rerun-$STAMP
    git tag -f $TAG
    say "freeze: $TAG = $(git rev-parse --short HEAD)"
    # clean rebuild of the main install from exactly this commit
    (cd src && wclean > /dev/null 2>&1 && nice -n 19 wmake -j 8 libso > /tmp/fr_lib.log 2>&1) || { say "lib build failed"; exit 1; }
    for a in applications/coupledFoam applications/test/Test-* applications/utilities/*; do
        [ -d "$a/Make" ] || continue
        (cd $a && wclean > /dev/null 2>&1; nice -n 19 wmake > /tmp/fr_app.log 2>&1) || { say "build failed: $a"; exit 1; }
    done
    say "freeze: build ok ($(ldd $(which coupledFoam) | grep -o '/home/jonas/OpenFOAM[^ ]*libcoupledFoam.so'))"
fi

if has archive; then
    A=$R/run/_archive_$STAMP
    mkdir -p $A
    # move coupledFoam run dirs; keep meshes (*_mesh*), simpleFoam references
    # (ref_*), unit_cavity and the archive dirs themselves
    for d in $R/run/*/; do
        n=$(basename $d)
        case $n in *_mesh*|ref_*|unit_cavity|_archive_*) continue;; esac
        mv "$d" "$A/"
    done
    RA=$R/results/_archive_$STAMP
    mkdir -p $RA
    # results/tests and results/gates are untracked run outputs; exploratory/ stays
    for sub in tests gates bench; do [ -d results/$sub ] && mv results/$sub $RA/$sub; done
    mkdir -p results/tests
    say "archive: run -> $A, results -> $RA"
fi

heavy_lane() {
    say "heavy lane start"
    for t in "test_T4_motorBike.py::test_T4[a]" "test_T4_motorBike.py::test_T4[b]" "test_T5_ahmed.py::test_T5"; do
        say "heavy: $t"
        CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 CF_T5_MESH=coarse nice -n 19 $PY -m pytest "tests/$t" --heavy -q -p no:cacheprovider \
            > $R/run/final_heavy_$(echo $t | tr -c 'A-Za-z0-9' '_').log 2>&1
        say "heavy: $t -> $(grep -E 'passed|failed' $R/run/final_heavy_$(echo $t | tr -c 'A-Za-z0-9' '_').log | tail -1)"
    done
    say "heavy lane done"
}

light_lane() {
    say "light lane start (CPUs $LIGHT_CPUS)"
    for t in test_unit.py test_env.py test_diagnostics.py test_T0_cavity.py test_T1_pitzDaily.py \
             test_T2_backwardFacingStep.py test_T3_airFoil.py test_restart.py test_fpe.py; do
        CF_MPI_BIND=none taskset -c $LIGHT_CPUS nice -n 19 $PY -m pytest tests/$t --ranks 1,4 -q -p no:cacheprovider -m "not heavy" \
            > $R/run/final_light_${t%.py}.log 2>&1
        say "light: $t -> $(grep -E 'passed|failed|no tests ran' $R/run/final_light_${t%.py}.log | tail -1)"
    done
    # small benchmark (np1, 3 repetitions) on the light cores while the heavy lane runs
    CF_MPI_BIND=none taskset -c $LIGHT_CPUS nice -n 19 $PY bench/run_bench.py --cases T1,T2,T3-SST,T3-GEKO --repeats 3 --allow-busy > $R/run/final_bench_small.log 2>&1
    say "light: benchmark small rc $?"
    say "light lane done"
}

if has heavy && has light; then
    heavy_lane & hp=$!
    light_lane & lp=$!
    wait $hp $lp
elif has heavy; then heavy_lane
elif has light; then light_lane
fi

if has scaling; then
    # ranks and case per user decision (see OPUS_TASKS / D-0xx); default T4b
    # runs alone on the machine: up to 16 physical cores (user budget 16/30)
    # D-059: T4a, 150 iterations, physical cores only
    export CF_SCALING_RANKS=${CF_SCALING_RANKS:-1,2,4,8,12,16} CF_SCALING_MESH=a CF_SCALING_ITERS=150
    say "scaling start (ranks $CF_SCALING_RANKS)"
    CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 CF_MPI_BIND=core nice -n 19 $PY -m pytest tests/test_scaling.py --heavy -q -p no:cacheprovider > $R/run/final_scaling.log 2>&1
    say "scaling -> $(grep -E 'passed|failed' $R/run/final_scaling.log | tail -1)"
fi

if has bench; then
    say "benchmark start"
    # D-059: heavy cases 1 repetition; T4a every configuration, T4b/T5 only E
    CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 nice -n 19 $PY bench/run_bench.py --cases T4a --no-scope --repeats 1 > $R/run/final_bench_T4a.log 2>&1
    say "benchmark T4a rc $?"
    CF_HEAVY_NP=10 CF_MPI_CPUSET=0-9 CF_T5_MESH=coarse nice -n 19 $PY bench/run_bench.py --cases T4b,T5 --configs E --repeats 1 > $R/run/final_bench_T4bT5.log 2>&1
    say "benchmark T4b/T5 rc $?"
fi

if has report; then
    SHA=$(git rev-parse HEAD)
    say "report (guard mode, commit $SHA)"
    nice -n 19 $PY bench/make_report.py --commit $SHA --allow-historical > $R/run/final_report.log 2>&1
    (cd report/paper && nice -n 19 make > $R/run/final_make.log 2>&1)
    say "report: $(grep -ci pending $R/run/final_report.log) pending notes; pdf rc $?"
fi
say "final re-run finished: $phases"
