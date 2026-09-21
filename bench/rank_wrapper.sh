#!/bin/sh
#------------------------------------------------------------------------------
# rank_wrapper.sh - run one MPI rank (or a serial run) under /usr/bin/time -v
# and keep the report per rank (spec 14): peak RSS, user and system CPU time,
# wall time.
#
# Used through the Allrun hook CF_RANK_WRAPPER:
#     CF_RANK_WRAPPER="/path/bench/rank_wrapper.sh" ./Allrun ...
# Reports go to $CF_TIMING_DIR (default ./timing) as
#     <app>.rank<N>.time   (GNU time -v output)
#------------------------------------------------------------------------------
dir="${CF_TIMING_DIR:-timing}"
mkdir -p "$dir"
rank="${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}"
app="$(basename "$1")"
exec /usr/bin/time -v -o "$dir/$app.rank$rank.time" "$@"
