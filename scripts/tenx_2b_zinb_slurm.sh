#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH -p LM
#SBATCH --mem 500GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=diyadas@berkeley.edu,rchance@berkeley.edu

module load gcc

ncores=$1
expt="regenK5"
#expt="ob"
#expt="cortex"
idfilt="FALSE"

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_2b_zinb_slurm_${job}_${NOW}.Rout"
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
       echo "usage: tenx_2b_zinb.sh ncores" >&2
       exit 2
}

[[ $# -eq 1 ]] || usage  # fail if incorrect number of args and print usage info

while true; do free -h >> 'memorylogs/tenx_2b_zinb_slurm_'$NOW'_memory.out'; sleep 15; done &

run tenx_2b_zinb.R \
   env R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --ncores "$ncores" --idfilt "$idfilt"

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2b_zinb.R --args --expt regen --ncores $ncores > 'tenx_2b_zinb'$NOW'.Rout' 2>&1
#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2b_zinb.R --args --expt ob  --ncores $ncores > 'tenx_2b_zinb'$NOW'.Rout' 2>&1
#write 24 and then 2
