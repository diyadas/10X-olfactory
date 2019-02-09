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
NOW=$(date +"_%Y%m%d-%H%M%S")

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --norm "none,fq,ruv_k=1,no_bio,no_batch" --method scone --clusmethod snn > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt regen --ncores $ncores --norm $2 --method scone > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1


#  R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt ob --ncores $ncores --norm $2 \
#  --method scone \
#  --markerfile "OBmarkers.txt" \
#  #'ob_scone_none,fq,ruv_k=3,no_bio,batch_res05_DE_OneAgainstAll_20181220_192548_top5abslogFC.txt' \
# > 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_viz_clusterPlots.R --args --expt ob --ncores $ncores --norm $2 \
#--method scone \
#--markerfile "OBmarkers.txt" \
#--clusmethod rsec \
#'ob_scone_none,fq,ruv_k=3,no_bio,batch_res05_DE_OneAgainstAll_20181220_192548_top5abslogFC.txt' \
#> 'tenx_viz_clusterPlots'$NOW'.Rout' 2>&1
 