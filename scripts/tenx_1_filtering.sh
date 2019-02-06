#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M diyadas@berkeley.edu
#$ -m beas
#

ncores=23
NOW=$(date +"_%Y%m%d-%H%M%S")

R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args \
--expt regen --ncores $ncores --aggr regen4 \
--exptinfo regen_exptinfo.csv \
--runQC TRUE \
--exclude double \
> 'tenx_1_filtering'$NOW'.Rout' 2>&1

#R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla < tenx_1_filtering.R --args \
#--expt ob --ncores $ncores --aggr RCOB2AB56 \
#--exptinfo rcob_exptinfo.csv \
#--hkfile hkpackage \
#--posctrlfile "../ref/OBmarkers.txt" \
#--runQC TRUE \
#--exclude aon \
#> 'tenx_1_filtering'$NOW'.Rout' 2>&1


####
  # make_option("--expt", type = "character", help = "Experiment ID"),
  # make_option("--ncores", default = "1", type = "double"),
  # make_option("--aggr", type = "character",
  #             help = "name of aggregate ID/directory"),
  # make_option("--exptinfo", type = "character",
  #             help = "text file with list of cellranger experiments, from
  #             aggregate...csv, expanded with columns for expt and batch: \n
  #             column 1: library_id, column 3: expt, column 4: batch"),
  # make_option("--annotation", default = "GRCm38p4Mm10", type = "character",
  #             help = "name of genomic reference used by cellranger"),
  # make_option("--hkfile", default = "../ref/hkl615.txt", type = "character",
  #             help = "path to file containing housekeeping genes"),
  # make_option("--posctrlfile", default = "../ref/oeRegPosCon.txt",
  #             type = "character",
  #             help = "path to file containing positive control genes"),
  # make_option("--runQC", default = FALSE, type = "logical",
  #             help = "whether to calculate QC metrics, will error if
  #             QC metrics have not previously been calculated"),
  # make_option("--fast", default = FALSE, type = "logical",
  #             help = "whether to use fast (approximate) PCA")
  # )

### 