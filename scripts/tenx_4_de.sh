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
method=$2 #zinb or scone
normalization=$3
expt="ob"
#expt="cortex"
idfilt="TRUE"

job=$(basename "$0")
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_4_de_${job}_${NOW}.Rout"
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
       echo "usage: tenx_4_de.sh ncores method normalization" >&2
       exit 2
}

[[ $# -eq 3 ]] || usage  # fail if incorrect number of args and print usage info


while true; do free -h >> 'tenx_4_de_'$NOW'_memory.out'; sleep 15; done &


run tenx_4_de.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --normalization "$normalization" \
       --ncores "$ncores" --method "$method" 


#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_4_de.R --args --expt ob --ncores $ncores --normalization "none,fq,ruv_k=3,no_bio,batch" --method scone > 'tenx_4_de'$NOW'.Rout' 2>&1
