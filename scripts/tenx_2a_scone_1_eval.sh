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

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_1_eval.R --args --expt regen --ncores $ncores --subsample TRUE > 'tenx_2a_scone_1_eval'$NOW'.Rout'

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_1_eval.R --args --expt RCOB  --ncores $ncores --subsample TRUE > 'tenx_2a_scone_1_eval'$NOW'.Rout'