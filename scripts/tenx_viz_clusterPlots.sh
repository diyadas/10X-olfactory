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
#expt="cortex"
#markerfile="cortexgenes.txt"
expt="ob"
markerfile="OBmarkers.txt"
#expt="regen"
#markerfile="oe_markers32+regen.txt"
idfilt="FALSE"

job=$(basename "$0")
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_viz_clusterPlots_${job}_${NOW}.Rout"
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
       echo "usage: tenx_viz_clusterPlots.sh ncores method normalization" >&2
       exit 2
}

[[ $# -eq 3 ]] || usage  # fail if incorrect number of args and print usage info


while true; do free -h >> 'tenx_viz_clusterPlots_'$NOW'_memory.out'; sleep 15; done &


run tenx_viz_clusterPlots.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --normalization "$normalization" \
       --ncores "$ncores" --method "$method" --markerfile "$markerfile"


#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --norm "none,fq,ruv_k=1,no_bio,no_batch" --method scone --clusmethod rsec > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --norm $2 --method scone > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1


#  R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt ob --ncores $ncores --norm $2 \
#  --method scone \
#  --markerfile "OBmarkers.txt" \
  #'ob_scone_none,fq,ruv_k=3,no_bio,batch_res05_DE_OneAgainstAll_20181220_192548_top5abslogFC.txt' \
# > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt ob --ncores $ncores --norm $2 \
#--method scone \
#--markerfile "OBmarkers.txt" \
#--clusmethod snn \
#--samplesort dendrogramValue \
#--seures "res.2" \
#> 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1
 