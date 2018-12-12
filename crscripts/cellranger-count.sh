#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas

NOW=$(date +'%Y%M%D_%H%M%S')
sample=$1

cellranger count --id=$sample --sample=$sample --fastqs=$SHARE/10X-olfactory/output/cellranger/HWVK5BBXX/outs/fastq_path/HWVK5BBXX --transcriptome=$SHARE/annotations/GRCm38p4Mm10