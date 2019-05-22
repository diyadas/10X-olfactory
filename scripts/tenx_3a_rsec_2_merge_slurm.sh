#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH -p LM
#SBATCH --mem 500GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=diyadas@berkeley.edu,rchance@berkeley.edu

 
ncores=$1
method=$2 #zinb or scone
normalization=$3
#expt="cortex"
#markerfile="cortexgenes.txt"
#expt="ob"
#markerfile="OBmarkers.txt"
#whichmerge= "plotdendro_adjP_cutoff_0.01"
seures="res.0.5"
#expt="regen"
expt="regenK5"
markerfile="oe_markers_regen.txt"
idfilt="FALSE"


job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_3a_rsec_2_merge_${job}_${NOW}.Rout"
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
       echo "usage: tenx_3a_rsec_2_merge_slurm.sh ncores method normalization" >&2
       exit 2
}

[[ $# -eq 3 ]] || usage  # fail if incorrect number of args and print usage info

while true; do free -h >> 'memorylogs/'$expt'_3a_rsec_2_merge_'$NOW'_memory.out'; sleep 15; done &

module load gcc

run tenx_3a_rsec_2_merge.R \
    env R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla --args \
        --expt "$expt" \
        --ncores $ncores \
        --normalization  "$normalization" \
        --method "$method"  \
        --markerfile "$markerfile" \
        --idfilt "$idfilt" \
        --whichmerge "$whichmerge"\
        --samplesort primaryCluster


