# Normalization of 10X data using zinbwave
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Oct 15 19:07:35 2018

# Load command-line arguments
rm(list = ls())
options(getClass.msg = FALSE)
library(optparse)
option_list <- list(
  make_option("--expt", type = "character", help = "Experiment ID"),
  make_option("--ncores", default = "1", type = "double")
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

load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
logfiltCounts <- log2(assay(se_filtered)+1)
vars <- rowVars(logfiltCounts)
names(vars) <- rownames(se_filtered)
vars <- sort(vars, decreasing = TRUE)
vargenes <- names(vars)[1:1000]
zinb <- zinbFit(se_filtered[vargenes,], X = "~ log10_total_counts + pct_counts_in_top_200_features + ribo_pct", K = 10, epsilon = 1000)
zinb_obj <- zinbwave(se_filtered[vargenes[1:1000],], fitted_model = zinb, K = 10, epsilon = 1000)
save(zinb_obj, 
       file = file.path(outdir, pasteu(exptstr, "zinb.Rda")))
