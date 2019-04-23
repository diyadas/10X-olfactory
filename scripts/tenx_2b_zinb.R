# Normalization of 10X data using zinbwave
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Oct 15 19:07:35 2018

# Load command-line arguments
rm(list = ls())
options(getClass.msg = FALSE)
library(optparse)
option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--idfilt", type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
print(paste("Using ", opt$ncores, "cores"))

# Set up packages and parallel environment
library(BiocParallel)
register(MulticoreParam(workers = opt$ncores))
library(zinbwave)
library(Rtsne)
# Source helper functions
source("tenx_helper.R")

if (opt$idfilt) {
  idfiltstr <- "idfiltyes"
} else {
  idfiltstr <- "idfiltno"
}

print(opt)
mytimestamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz="America/Los_Angeles")
print(paste("Files produced by this script will be timestamped:", mytimestamp))

datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "1_se_filtqc", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

logfiltCounts <- log2(assay(se_filtered)+1)
vars <- rowVars(logfiltCounts)
names(vars) <- rownames(se_filtered)
vars <- sort(vars, decreasing = TRUE)

zinb <- zinbFit(se_filtered[names(vars)[1:1000],], X = "~ log10_total_counts + pct_counts_in_top_500_features + ribo_pct", K = 20, epsilon = 1000)
zinbparams <- "zinbFit(se_filtered[names(vars)[1:1000],], X = '~ log10_total_counts + pct_counts_in_top_500_features + ribo_pct', K = 20, epsilon = 1000)"
save(zinb, zinbparams, file = file.path(datdir, pasteu(exptstr,"zinb_data.Rda")))

W <- getW(zinb)
d <- dist(W)

save(W, zinbparams, file = file.path(datdir, pasteu0(exptstr, "zinbW", format(Sys.time(), "%Y%m%d_%H%M%S"), ".Rda")))

zinb_obj <- zinbwave(se_filtered[names(vars)[1:1000],], fitted_model = zinb, K = 20, epsilon = 1000)
save(zinb_obj, zinbparams, zinb, 
       file = file.path(datdir, pasteu(exptstr,idfiltstr, "zinb_data.Rda")))
