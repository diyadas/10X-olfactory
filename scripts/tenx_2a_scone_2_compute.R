# Normalization of 10X data using SCONE
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Fri Mar  1 10:57:31 2019

# Load command-line arguments
rm(list = ls())
options(getClass.msg = FALSE)
library(optparse)
option_list <- list(
  make_option("--expt", type = "character", help = "Experiment ID"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--normalization", type = "character"),
  make_option("--idfilt", type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")

print(opt$expt)
print(opt$ncores)
print(opt$normalization)

# Set up packages and parallel environment
library(BiocParallel)
register(MulticoreParam(workers = opt$ncores))
library(scone)

# Source helper functions
source("tenx_helper.R")

scone_fx <- function(SconeExperimentObj, ...){
  scone(SconeExperimentObj,
        scaling = list(none = identity,sum = SUM_FN, tmm = TMM_FN, uq = UQ_FN,
                       fq = FQT_FN, deseq = DESEQ_FN),
        k_ruv = 3,
        k_qc = 5,
        adjust_bio = "no",
        adjust_batch = "yes",
        eval_kclust = 5:15,
        zero = "postadjust", ... 
  )
}

load(file = file.path(datdir, pasteu0(exptstr, "2a_scone_1_eval_select", idfiltstr, ".Rda")))
scone_obj@scone_params <- scone_obj@scone_params[opt$normalization,]
scone_obj <- scone_fx(scone_obj, run = TRUE, return_norm = "in_memory")
save(scone_obj, file = file.path(datdir, pasteu0(exptstr, "scone", opt$normalization, "data", idfiltstr, ".Rda")))
