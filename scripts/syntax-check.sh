#!/bin/bash
#------------------------------------------------------------------------------
# syntax-check.sh - compile translation units with -fsyntax-only
#
# Catches typos and signature errors without a full wmake build. Uses the
# same warning flags as the project Make/options (-Werror included) and the
# OpenFOAM headers of the currently sourced environment as -isystem.
#
# Usage:
#   scripts/syntax-check.sh [-p DP|SPDP|SP] file.C [file.C ...]
#
# -p overrides -DWM_<precision>, e.g. "-p SP" checks that the code also
# compiles with scalar == float (the block library must be precision-agnostic,
# DECISIONS.md D-001).
#
# Requires an OpenFOAM environment (WM_PROJECT_DIR set). Runs niced, one file
# at a time, so it is safe on a busy machine.
#------------------------------------------------------------------------------
set -e

if [ -z "$WM_PROJECT_DIR" ]
then
    echo "No OpenFOAM environment sourced" >&2
    exit 2
fi

precision="$WM_PRECISION_OPTION"
if [ "$1" = "-p" ]
then
    precision="$2"
    shift 2
fi

root="$(cd "$(dirname "$0")/.." && pwd)"
src="$WM_PROJECT_DIR/src"

sysInc=(
    "$src/OpenFOAM/lnInclude"
    "$src/OSspecific/${WM_OSTYPE:-POSIX}/lnInclude"
    "$src/finiteVolume/lnInclude"
    "$src/meshTools/lnInclude"
    "$src/sampling/lnInclude"
    "$src/TurbulenceModels/turbulenceModels/lnInclude"
    "$src/TurbulenceModels/incompressible/lnInclude"
    "$src/transportModels"
    "$src/transportModels/incompressible/lnInclude"
    "$src/transportModels/incompressible/singlePhaseTransportModel"
    "$src/functionObjects/forces/lnInclude"
)

projInc=(
    "$root/src/blockMatrix"
    "$root/src/blockSolvers"
    "$root/src/blockSolvers/blockSmoothers"
    "$root/src/blockSolvers/blockPreconditioners"
    "$root/src/assembly"
    "$root/src/control"
    "$root/src/io"
    "$root/src/include"
)

flags=(
    -std=c++17 -m64 -pthread
    -DOPENFOAM=2606 "-DWM_${precision}" "-DWM_LABEL_SIZE=${WM_LABEL_SIZE:-32}"
    -DNoRepository -ftemplate-depth=1000
    -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor
    -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes
    -Wno-unknown-pragmas
    -Wfloat-conversion -Werror
    -fno-math-errno
    -fsyntax-only
)

for d in "${sysInc[@]}"; do flags+=(-isystem "$d"); done
for d in "${projInc[@]}"; do flags+=(-I "$d"); done

status=0
for f in "$@"
do
    dir="$(cd "$(dirname "$f")" && pwd)"
    if nice -n 19 g++ "${flags[@]}" -I "$dir" "$f"
    then
        echo "OK    $f"
    else
        echo "FAIL  $f"
        status=1
    fi
done

exit $status

#------------------------------------------------------------------------------
