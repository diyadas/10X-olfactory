#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas
#

ncores=$1
normalization=$2
expt="ob"
#expt="cortex"
idfilt="FALSE"

job=$(basename "$0")
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_2a_scone_2_compute_${job}_${NOW}.Rout"
exec >> "$LOG" 2>&1 || exit 1     # redirect stdout and error to log file, will fail if the logs directory doesn't exist

# as refactored by JW Adams from last commit
# usage: run input_file command to run
run() {
       local input="$1"
       shift;  # remove file from argument list
       echo "command: $* < $input"
       echo "--------------"
       "$@" < "$input"
}

usage() {
       echo "usage: tenx_2a_scone_2_compute.sh ncores normalization" >&2
       exit 2
}

[[ $# -eq 2 ]] || usage  # fail if incorrect number of args and print usage info

while true; do free -h >> 'tenx_2b_scone_2_compute_'$NOW'_memory.out'; sleep 15; done &

run tenx_2a_scone_2_compute.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --ncores "$ncores" --idfilt "$idfilt"  --normalization "$normalization"

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_2_compute.R --args --expt regen --ncores $ncores --normalization $2  \
#> 'tenx_2a_scone_2_compute'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_2_compute.R --args --expt ob  --ncores $ncores --normalization $2  \
#> 'tenx_2a_scone_2_compute'$NOW'.Rout' 2>&1