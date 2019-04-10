#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu,rchance@berkeley.edu
#$ -m beas

NOW=$(date +'%Y%M%D_%H%M%S')

#cellranger count --id=$sample --sample=$sample --fastqs=$SHARE/10X-olfactory/output/cellranger/HWVK5BBXX/outs/fastq_path/HWVK5BBXX --transcriptome=$SHARE/annotations/GRCm38p4M10
#FASTQPATH="/share/groups/diya-russell/10X-olfactory/output/seq/Stafford/RBP4B"
#FASTQPATH="/share/groups/diya-russell/10X-olfactory/output/seq/190315_100PE_HS4K2B"
#FASTQPATH="/share/groups/diya-russell/10X-olfactory/output/seq/190322_100PE_HS4K2A" #wholeOE mini injury timecourse
FASTQPATH="/share/groups/diya-russell/10X-olfactory/output/seq/regen_whOE_batch2" #wholeOE mini injury timecourse batch2
cellranger count --id=$1 --sample=$2 --fastqs=$FASTQPATH --transcriptome=$SHARE/annotations/GRCm38p4Mm10

#id regen_24HPI-2 this is a folder -1 is RF  -2 is RKC   sample SI-GA-G4-a  