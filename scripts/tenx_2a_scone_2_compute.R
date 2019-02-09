# Normalization of 10X data using SCONE
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Wed Jul  4 23:52:48 2018

# Load command-line arguments
rm(list = ls())
options(getClass.msg = FALSE)
library(optparse)
option_list <- list(
  make_option("--expt", type = "character", help = "Experiment ID"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--normalization", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
outdir <- file.path("../output", exptstr, "data")
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

load(file = file.path(outdir, pasteu(exptstr, "scone_eval_select.Rda")))
scone_obj@scone_params <- scone_obj@scone_params[opt$normalization,]
scone_obj <- scone_fx(scone_obj, run = TRUE, return_norm = "in_memory")
save(scone_obj, file = file.path(outdir, pasteu(exptstr, "scone", opt$normalization, "data.Rda")))
