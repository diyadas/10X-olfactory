#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas
#

NOW=$(date +"_%m%d%Y-%H%M%S")

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt regen --ncores $1 --normalization "none,fq,ruv_k=2,no_bio,no_batch" --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt regen --ncores $1 --normalization $2 --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt ob --ncores $1 --normalization $2 --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_3b_snn.R --args --expt ob --ncores $1 --normalization "none,fq,ruv_k=3,no_bio,batch"  --method scone > 'tenx_3b_snn'$NOW'.Rout' 2>&1