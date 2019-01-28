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

while true; do free -h >> 'tenx_3_rsec'$NOW'_memory.out'; sleep 15; done &

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3a_rsec.R --args --expt ob --ncores $ncores --normalization $2 --method scone > 'tenx_3a_rsec'$NOW'.Rout' 2>&1 
