#!/bin/bash
# Build coupledFoam (library, solver, test apps, utilities) of the current
# checkout into the SP user install. Run inside the SP environment:
#   /home/jonas/bin/cfenv sp bash bench/build_sp_install.sh
# The SP user install lives under the WM_OPTIONS linux64GccSPInt32Opt
# platform directory, separate from the DP one (D-062).
set -u
R=${CF_REPO:-/home/jonas/coupledFoam}
[ "$WM_PRECISION_OPTION" = "SP" ] || { echo "not an SP environment ($WM_OPTIONS)"; exit 2; }
cd $R/src && wclean > /dev/null 2>&1
nice -n 19 wmake -j ${WM_NCOMPPROCS:-8} libso > /tmp/sp_lib.log 2>&1 || {
    echo "SP LIB FAILED"; grep -m12 -iE "error" /tmp/sp_lib.log; exit 1; }
for a in $R/applications/coupledFoam $R/applications/test/Test-* $R/applications/utilities/*; do
    [ -d "$a/Make" ] || continue
    (cd "$a" && wclean > /dev/null 2>&1; nice -n 19 wmake > /tmp/sp_app.log 2>&1) || {
        echo "SP APP FAILED: $a"; grep -m10 -iE "error" /tmp/sp_app.log; exit 1; }
done
echo "SP build ok: $(ldd $(which coupledFoam) | grep -o '/home/jonas/OpenFOAM[^ ]*libcoupledFoam.so')"
Test-precision 2>&1 | tail -3
