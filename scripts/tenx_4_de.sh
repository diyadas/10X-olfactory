#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas
#

ncores=$23
NOW=$(date +"_%Y%m%d-%H%M%S")

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_4_de.R --args --expt ob --ncores $ncores --normalization "none,fq,ruv_k=3,no_bio,batch" --method scone > 'tenx_4_de'$NOW'.Rout' 2>&1
