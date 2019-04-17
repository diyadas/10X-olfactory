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

#cellranger mkfastq --run=$SHARE/10X-olfactory/output/ob/seq/181015_K00364_0106_BHWVK5BBXX --csv=$SHARE/10X-olfactory/output/ob/RCOB6-samplesheet.csv --ignore-dual-index --delete-undetermined
#cellranger mkfastq --run=$SHARE/10X-olfactory/output/seq/190322_100PE_HS4K2A --csv=$SHARE/10X-olfactory/output/seq/190322_100PE_HS4K2A/samplesheet_190322.csv --ignore-dual-index --delete-undetermined
cellranger mkfastq --run=$SHARE/10X-olfactory/output/seq/regen_whOE_batch2 --csv=$SHARE/10X-olfactory/output/seq/regen_whOE_batch2/samplesheet_whOE_batch2.csv --ignore-dual-index --delete-undetermined
# --id=$sample --fastqs=$SHARE/10X-olfactory/output/ob/fastq/RCOB6 --transcriptome=$SHARE/annotations/GRCm38p4Mm10