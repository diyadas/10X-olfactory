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
expt="ob"
normalization=$2

NOW=$(date +"_%Y%m%d-%H%M%S")

CMD='R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_2_compute.R --args --expt '$expt'  --ncores '$ncores' --normalization '$normalization

echo "The command invoked: " >> logs/$expt"_tenx_2a_scone_2_compute"$NOW".Rout" 2>&1
echo $CMD >> logs/$expt"_tenx_2a_scone_2_compute"$NOW".Rout" 2>&1
echo "--------------" >> logs/$expt"_tenx_2a_scone_2_compute"$NOW".Rout" 2>&1
eval $CMD >> logs/$expt"_tenx_2a_scone_2_compute"$NOW".Rout" 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_2_compute.R --args --expt regen --ncores $ncores --normalization $2  \
#> 'tenx_2a_scone_2_compute'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_2a_scone_2_compute.R --args --expt ob  --ncores $ncores --normalization $2  \
#> 'tenx_2a_scone_2_compute'$NOW'.Rout' 2>&1