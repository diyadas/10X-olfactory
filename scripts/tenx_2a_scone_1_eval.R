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
  make_option("--subsample", default = TRUE, type = "logical")
)
opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
outdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")

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

load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
scone_obj <- SconeExperiment(se_filtered,
                             which_qc = 
                               which(colnames(colData(se_filtered)) %in% 
                                       colnames(qc)),
                             which_negconruv = 5L,
                             which_negconeval = 4L,
                             which_poscon = 3L,
                             which_batch = 1L,
                             which_bio = 2L
)
scone_params <- scone_fx(scone_obj, run = FALSE)

if (opt$subsample) {
  print("Running scone with subsampling")
  N <- NCOL(scone_obj)
  n <- ceiling(N / 10)
  sub_ranks <- replicate(10, {
    print("New replicate...")
    idx <- sample(N, n)
    scone_sub <- scone_obj[, idx]
    scone_sub <- scone_sub[rowSums(assay(scone_sub)) > 0, ]     # filter genes
    scone_sub <- scone_fx(scone_sub)     # run scone
    return(list(get_scores(scone_sub), get_score_ranks(scone_sub)))
  } , simplify = FALSE)
  save(scone_obj, sub_ranks, 
       file = file.path(outdir, pasteu(exptstr, "scone_eval_full.Rda")))
  
  scone_ranks <- sort(rowMeans(
    sapply(sub_ranks, function(x) x[[2]][names(sub_ranks[[1]][[2]])])), 
    decreasing = TRUE)
  scone_scores <- lapply(sub_ranks, function(x) x[[1]][rownames(sub_ranks[[1]][[1]]),])
  scone_scores_mean <- Reduce("+",scone_scores)/length(scone_scores)
  scone_scores_mean <- scone_scores_mean[names(scone_ranks),]  
  scone_obj@scone_params <- scone_params@scone_params[names(scone_ranks)[1:8],] 
  save(scone_obj, scone_ranks, scone_scores, scone_scores_mean,
       file = file.path(outdir, pasteu(exptstr, "scone_eval_select.Rda")))
} else {
  print("Running scone without subsampling")
  scone_obj <- scone_fx(scone_obj)
  save(scone_obj, 
       file = file.path(outdir, pasteu(exptstr, "scone_eval_full.Rda")))
  scone_obj@scone_params <- scone_obj@scone_params[1:6,]
  save(scone_obj, 
       file = file.path(outdir, pasteu(exptstr, "scone_eval_select.Rda")))
}