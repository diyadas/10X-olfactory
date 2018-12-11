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
print(paste("Using ", opt$ncores, "cores"))

# Set up packages and parallel environment
library(BiocParallel)
register(MulticoreParam(workers = opt$ncores))
library(zinbwave)
library(Rtsne)
# Source helper functions
source("tenx_helper.R")

load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
logfiltCounts <- log2(assay(se_filtered)+1)
vars <- rowVars(logfiltCounts)
names(vars) <- rownames(se_filtered)
vars <- sort(vars, decreasing = TRUE)

zinb <- zinbFit(se_filtered[names(vars)[1:1000],], X = "~ log10_total_counts + pct_counts_in_top_200_features + ribo_pct", K = 20, epsilon = 1000)
zinbparams <- "zinbFit(se_filtered[names(vars)[1:1000],], X = '~ log10_total_counts + pct_counts_in_top_200_features + ribo_pct', K = 20, epsilon = 1000)"
save(zinb, zinbparams, file = file.path(outdir, paste0(exptstr,"_zinb.Rda")))

W <- getW(zinb)
d <- dist(W)

save(W, zinbparams, file = file.path(outdir, pasteu(exptstr, "zinbW", Sys.Date(), ".Rda")))

zinb_obj <- zinbwave(se_filtered[names(vars)[1:1000],], fitted_model = zinb, K = 20, epsilon = 1000)
save(zinb_obj, zinbparams,
       file = file.path(outdir, pasteu(exptstr, "zinb.Rda")))
