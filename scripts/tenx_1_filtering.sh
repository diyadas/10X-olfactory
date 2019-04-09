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
#expt="cortex"
#aggr="Rbp4AB"
#exptinfo="cortex_exptinfo.csv"
#posctrlfile="../ref/cortexgenes.txt"


expt="ob"
aggr="RCOB2AB567"
exptinfo="rcob_exptinfo.csv"
posctrlfile="../ref/OBmarkers.txt"

hkfile="hkpackage"

runQC="TRUE"
exclude=""

job=$(basename "$0")
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="logs/${expt}_1_filtering_${job}_${NOW}.Rout"
exec >> "$LOG" 2>&1 || exit 1     # redirect stdout and error to log file, will fail if the logs directory doesn't exist

# as refactored by JW Adams from last commit
# usage: run input_file command to run
run() {
       local input="$1"
       shift;  # remove file from argument list
       echo "command: $* < $input"
       echo "--------------"
       "$@" < "$input"
}

usage() {
       echo "usage: tenx_1_filtering.sh ncores" >&2
       exit 2
}

[[ $# -eq 1 ]] || usage  # fail if incorrect number of args and print usage info

run tenx_1_filtering.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla --args \
       --expt "$expt" --ncores "$ncores" --aggr "$aggr" \
       --exptinfo "$exptinfo" \
       --hkfile "$hkfile" \
       --posctrlfile "$posctrlfile" \
       --runQC "$runQC" \
       --exclude "$exclude"


#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args \
#--expt $EXPT --ncores $ncores --aggr reg_whOE \
#--exptinfo regen_whoe_exptinfo.csv \
#--runQC TRUE \
#--exclude double \
#> logs/$EXPT'_tenx_1_filtering'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args \
#--expt $EXPT --ncores $ncores --aggr RCOB2AB567 \
#--exptinfo rcob_exptinfo.csv \
#--hkfile hkpackage \
#--posctrlfile "../ref/OBmarkers.txt" \
#--runQC TRUE \
#--exclude aon \
#> logs/$EXPT'_tenx_1_filtering'$NOW'.Rout' 2>&1

# R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args \
# --expt $EXPT  --ncores $ncores --aggr Rbp4AB \
# --exptinfo cortex_exptinfo.csv \
# --hkfile hkpackage \
# --posctrlfile "../ref/cortexgenes.txt" \
# --runQC TRUE \
# #--exclude gad1 \
# > logs/$EXPT'_tenx_1_filtering'$NOW'.Rout' 2>&1



# option_list <- list(
#   make_option("--expt", type = "character", help = "Experiment ID"),
#   make_option("--ncores", default = "1", type = "double"),
#   make_option("--aggr", type = "character",
#               help = "name of aggregate ID/directory"),
#   make_option("--exptinfo", type = "character",
#               help = "text file with list of cellranger experiments, from 
#               aggregate...csv, expanded with columns for expt and batch: \n 
#               column 1: library_id, column 3: expt, column 4: batch"),
#   make_option("--annotation", default = "GRCm38p4Mm10", type = "character",
#               help = "name of genomic reference used by cellranger"),
#   make_option("--hkfile", default = "../ref/hkl615.txt", type = "character", 
#               help = "path to file containing housekeeping genes"),
#   make_option("--posctrlfile", default = "../ref/oeRegPosCon.txt", 
#               type = "character",
#               help = "path to file containing positive control genes"),
#   make_option("--runQC", default = FALSE, type = "logical",
#             help = "whether to calculate QC metrics, will error if
# 	            QC metrics have not previously been calculated"),
#   make_option("--fast", default = FALSE, type = "logical",
#             help = "whether to use fast (approximate) PCA"),
#   make_option("--exclude", default = NULL, type="character", help = "name for excluded samples, if given")
#   )