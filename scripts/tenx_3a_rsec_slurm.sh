#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH -p LM
#SBATCH --mem 500GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=diyadas@berkeley.edu,rchance@berkeley.edu


#### #SBATCH -t 48:00:00 #SBATCH --mem 500GB

module load gcc
 
ncores=$1
method=$2 #zinb or scone
normalization=$3
expt="ob"
#expt="cortex"
idfilt="FALSE"

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_3a_rsec_${job}_${NOW}.Rout"
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
       echo "usage: tenx_3a_rsec.sh ncores method normalization" >&2
       exit 2
}

[[ $# -eq 3 ]] || usage  # fail if incorrect number of args and print usage info


while true; do free -h >> 'tenx_3a_rsec_'$NOW'_memory.out'; sleep 15; done &


run tenx_3a_rsec.R \
   env R_LIBS=/pylon5/ib5phhp/diyadas/rpack/3.5/  R --vanilla --args \
       --expt "$expt" --normalization "$normalization" \
       --ncores "$ncores" --method "$method" --idfilt "$idfilt" \
       --sequential TRUE \
       --subsample TRUE
