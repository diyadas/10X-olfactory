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

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2b_zinb.R --args --expt regen --ncores $ncores > 'tenx_2b_zinb'$NOW'.Rout' 2>&1
