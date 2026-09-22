#!/bin/bash
#------------------------------------------------------------------------------
# TASK 5 (D-045) acceptance, still to run in a quiet window:
#   (a) tests/test_diagnostics.py, np1 and np4
#   (b) unit battery, non-heavy: tests/test_unit.py (np1 + np4) and test_env
#   (c) T0 Re100 np1 wall-time medians: the pre-TASK-5 build (93dc500) vs.
#       this build at diagnostics level 0/1/2, 3 interleaved rounds
#       (gates: level 0 within +-5 % of the pre-TASK-5 build, level 1 <= 2 %,
#       level 2 <= 15 %)
#
# Usage (from any shell; sources the system OpenFOAM through cfenv):
#   /home/jonas/bin/cfenv sys bash bench/task5_acceptance.sh [a] [b] [c]
# With no argument all three run. Output: run/task5_acceptance.log and
# run/overhead.txt of this worktree. Needs ~4 free cores for (a)/(b).
#------------------------------------------------------------------------------
set -u
WT=/home/jonas/cf_diag                       # task5-diagnostics worktree
DIAG=/home/jonas/cf_diag_platform            # its private install
BASE=/home/jonas/cf_diag_base_platform       # pre-TASK-5 install (93dc500)
BASEWT=/home/jonas/cf_diag_base
PY=/home/jonas/OF/venv/bin/python
LOG=$WT/run/task5_acceptance.log
steps=${*:-a b c}
mkdir -p $WT/run
ORIGPATH=$PATH; ORIGLD=${LD_LIBRARY_PATH:-}

use() {  # use <install>
    export FOAM_USER_LIBBIN=$1/lib FOAM_USER_APPBIN=$1/bin
    export PATH=$1/bin:$ORIGPATH LD_LIBRARY_PATH=$1/lib:$ORIGLD
}

use $DIAG
echo "== $(date) steps: $steps, load $(cut -d' ' -f1-3 /proc/loadavg)" | tee -a $LOG
ldd "$(which coupledFoam)" | grep coupledFoam | tee -a $LOG
cd $WT || exit 1

case " $steps " in *" a "*)
    echo "== (a) test_diagnostics np1+np4" | tee -a $LOG
    CF_MPI_BIND=none $PY -m pytest tests/test_diagnostics.py --ranks 1,4 -q 2>&1 | tail -5 | tee -a $LOG
esac

case " $steps " in *" b "*)
    echo "== (b) unit battery (non-heavy) + test_env" | tee -a $LOG
    CF_MPI_BIND=none $PY -m pytest tests/test_unit.py tests/test_env.py -m "not heavy" -q 2>&1 | tail -8 | tee -a $LOG
esac

case " $steps " in *" c "*)
    echo "== (c) T0 Re100 np1 overhead" | tee -a $LOG
    if [ ! -x $BASE/bin/coupledFoam ]; then
        # pre-TASK-5 baseline build
        [ -d $BASEWT ] || git -C /home/jonas/coupledFoam worktree add --detach $BASEWT 93dc500
        ( use $BASE; mkdir -p $BASE/lib $BASE/bin
          cd $BASEWT/src && nice -n 19 wmake -j 2 libso > $BASE/build.log 2>&1
          cd $BASEWT/applications/coupledFoam && nice -n 19 wmake -j 2 >> $BASE/build.log 2>&1 )
    fi
    OUT=$WT/run/overhead.txt
    : > $OUT
    for r in 1 2 3; do
      for c in base l0 l1 l2; do
        if [ $c = base ]; then use $BASE; else use $DIAG; fi
        run=$WT/run/acc_${c}_$r
        rm -rf $run; cp -r $WT/cases/T0_cavity $run; cd $run
        case $c in l1|l2)
            foamDictionary -disableFunctionEntries -entry coupled/diagnostics \
                -set "{ level ${c#l}; }" system/fvSolution > /dev/null 2>&1 ;;
        esac
        CF_MPI_BIND=none nice -n 19 ./Allrun -solver coupledFoam -np 1 > log.Allrun 2>&1
        echo "$c $r rc=$? $(grep '^coupledFoam: iterations' log.coupledFoam | cut -d, -f1) $(grep '^coupledFoam: wall time' log.coupledFoam | cut -d'(' -f1) load=$(cut -d' ' -f1 /proc/loadavg)" | tee -a $OUT
      done
    done
    cd $WT
    $PY - "$OUT" <<'EOF' | tee -a $LOG
import re, statistics, sys
t = {}
for line in open(sys.argv[1]):
    c = line.split()[0]
    m = re.search(r"wall time ([\d.]+) s", line)
    if m:
        t.setdefault(c, []).append(float(m.group(1)))
med = {c: statistics.median(v) for c, v in t.items()}
b, l0 = med["base"], med["l0"]
print("medians [s]:", {c: round(v, 3) for c, v in med.items()})
print("level 0 vs base %+.2f %% (gate +-5)" % (100 * (l0 / b - 1)))
print("level 1 vs level 0 %+.2f %% (gate <= 2)" % (100 * (med["l1"] / l0 - 1)))
print("level 2 vs level 0 %+.2f %% (gate <= 15)" % (100 * (med["l2"] / l0 - 1)))
EOF
esac
echo "== done $(date)" | tee -a $LOG
