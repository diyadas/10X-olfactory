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
NOW=$(date +"_%m%d%Y-%H%M%S")

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --normalization $2 --method scone > 'tenx_viz_clusterPlots'$NOW'.Rout'
