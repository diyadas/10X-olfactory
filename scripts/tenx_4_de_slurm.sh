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
method=$2 #zinb or scone
normalization=$3
clusmethod=$4 #rsec or snn, if rsec add whichmerge

#expt="ob"
#expt="cortex"
#whichmerge="locfdr_cutoff_0.08"

expt="regen"
#expt="regenK5"
#whichmerge="adjP_cutoff_0.01"
#idfilt="TRUE"
idfilt="FALSE"

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_4_de_slurm_${job}_${NOW}.Rout"
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
       echo "usage: tenx_4_de.sh ncores method normalization clusmethod whichmerge" >&2
       exit 2
}

[[ $# -eq 4 ]] || usage  # fail if incorrect number of args and print usage info


while true; do free -h >> 'memorylogs/${expt}_tenx_4_de_slurm_'$NOW'_memory.out'; sleep 15; done &


run tenx_4_de.R \
   env R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --normalization "$normalization" --whichmerge "$whichmerge" \
       --ncores "$ncores" --method "$method" --clusmethod "$clusmethod"  --idfilt "$idfilt" 


#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_4_de.R --args --expt ob --ncores $ncores --normalization "none,fq,ruv_k=3,no_bio,batch" --method scone > 'tenx_4_de'$NOW'.Rout' 2>&1
