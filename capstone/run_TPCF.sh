#!/bin/bash
if [[ $# -ne 1 ]] ; then
    echo "usage: ./prog num_threads"
    exit 1
fi

export OMP_NUM_THREADS=$1
threadnums=$1

#time ./DM_CF_Parallel $1

'/usr/bin/time' -o time_stats_TPCF_omp${threadnums}.txt -a -f '%e %U %S' ./DM_CF_Parallel $threadnums 
