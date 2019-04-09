#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas

NOW=$(date +'%Y%M%D_%H%M%S')
#sample=$1

#cellranger count --id=$sample --sample=$sample --fastqs=$SHARE/10X-olfactory/output/cellranger/HWVK5BBXX/outs/fastq_path/HWVK5BBXX --transcriptome=$SHARE/annotations/GRCm38p4Mm10
cellranger aggr --id=regen_inclwhOE_20190407 --csv=$SHARE/10X-olfactory/output/regen/regen_exptinfo.csv --normalize=none
#cellranger aggr --id=RCOB2AB567 --csv=$SHARE/10X-olfactory/output/ob/rcob_exptinfo.csv --normalize=none
#cellranger aggr --id=cortex_Rbp4AB --csv=$SHARE/10X-olfactory/output/cortex/cortex_exptinfo.csv  --normalize=none
#cellranger aggr --id=regen_whOE_20190327 --csv=$SHARE/10X-olfactory/output/regen/crcount/reg_whOE/regen_whoe_exptinfo.csv  --normalize=none