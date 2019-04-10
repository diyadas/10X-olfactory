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
expt="ob"
method=$2
normalization=$3

job=$(basename "$0")
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_3b_snn_${job}_${NOW}.Rout"
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
       echo "usage: tenx_3b_snn.sh ncores method normalization" >&2
       exit 2
}

[[ $# -eq 3 ]] || usage  # fail if incorrect number of args and print usage info

run tenx_3b_snn.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --ncores "$ncores" --method "$method" --normalization "$normalization" 

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt regen --ncores $1 --normalization "none,fq,ruv_k=2,no_bio,no_batch" --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt regen --ncores $1 --normalization $2 --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt ob --ncores $1 --normalization $2 --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt ob --ncores $1 --normalization "none,fq,ruv_k=3,no_bio,batch"  --method zinb > 'tenx_3b_snn'$NOW'.Rout' 2>&1