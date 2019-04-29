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
normalization="$2"
#expt="cortex"
#markerfile="cortexgenes.txt"
#expt="ob"
#markerfile="OBmarkers.txt"
expt="regen"
markerfile="oe_markers32+regen.txt"
idfilt="FALSE"

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_edaviz_slurm_${job}_${NOW}.Rout"
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
       echo "usage: tenx_edaviz_slurm.sh ncores normalization" >&2
       exit 2
}

[[ $# -eq 2 ]] || usage  # fail if incorrect number of args and print usage info

run tenx_edaviz.R \
   env R_LIBS=/pylon5/ib5phhp/diyadas/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --ncores "$ncores" --normalization "$normalization" \
       --markerfile "$markerfile" --idfilt "$idfilt" 


#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_1_eval.R --args --expt regen --ncores $ncores --subsample TRUE > 'tenx_2a_scone_1_eval'$NOW'.Rout' 2>&1
#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_1_eval.R --args --expt ob  --ncores $ncores --subsample TRUE > 'tenx_2a_scone_1_eval'$NOW'.Rout' 2>&1
