#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH -p RM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=diyadas@berkeley.edu,rchance@berkeley.edu

module load gcc

ncores=$1
method=$2 #zinb or scone
normalization=$3
clusmethod=$4 #snn or rsec
#IF RSEC
merge=$5 #merged or notmerged
whichmerge=$6 #adjP_cutoff  or locfdr_cutoff

#expt="cortex"
#markerfile="cortexgenes.txt"
#expt="ob"
#markerfile="OBmarkers.txt"
seures="res.0.5"
expt="regen"
#expt="regenK5"
markerfile="oe_markers_regen.txt"
#markerfile="oe_markers32+regen.txt"
idfilt="FALSE"

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_viz_clusterPlots_slurm_${job}_${NOW}.Rout"
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
       echo "usage: tenx_viz_slurm_clusterPlots.sh ncores method normalization clusmethod" >&2
       exit 2
}

[[ $# -eq 4 ]] || usage  # fail if incorrect number of args and print usage info

while true; do free -h >> 'memorylogs/${expt}_tenx_viz_clusterPlots_slurm_'$NOW'_memory.out'; sleep 15; done &

#module load hdf5

#R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --norm "none,fq,ruv_k=1,no_bio,no_batch" --method scone --clusmethod snn > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

#R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla < \
#tenx_viz_clusterPlots.R --args \
#--expt regen \
#--ncores $ncores \
#--norm "none,fq,ruv_k=1,no_bio,no_batch" \
#--method scone \
#--clusmethod rsec \
#--samplesort primaryCluster > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

# tenx_viz_quicktsne.R --args \

run tenx_viz_clusterPlots.R \
  env R_LIBS=/pylon5/ib5phhp/shared/rpack/3.5/ R --vanilla --args \
--expt "$expt" \
--ncores "$ncores" \
--normalization  "$normalization" \
--method "$method"  \
--clusmethod "$clusmethod" \
--markerfile "$markerfile" \
--seures "$seures" \
--idfilt "$idfilt" \
--samplesort primaryCluster  
