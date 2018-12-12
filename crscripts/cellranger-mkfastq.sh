#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas

NOW=$(date +'%Y%M%D_%H%M%S')
module load bcl2fastq2-2.20

cellranger mkfastq --run=$SHARE/10X-olfactory/output/ob/seq/181015_K00364_0106_BHWVK5BBXX --csv=$SHARE/10X-olfactory/output/ob/RCOB6-samplesheet.csv --ignore-dual-index --delete-undetermined

# --id=$sample --fastqs=$SHARE/10X-olfactory/output/ob/fastq/RCOB6 --transcriptome=$SHARE/annotations/GRCm38p4Mm10