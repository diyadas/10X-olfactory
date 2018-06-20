#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu
#$ -m beas
#

ncores=19
NOW=$(date +"_%m%d%Y-%H%M%S")

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args --expt regen --ncores $ncores --aggr regen4 --exptinfo regen_exptinfo.csv > 'tenx_1_filtering'$NOW'.Rout'
