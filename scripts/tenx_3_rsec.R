# Clustering of 10X data using RSEC
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Jul  9 16:07:04 2018

# Load command-line arguments
library(scone)
library(clusterExperiment)
library(BiocParallel)
library(optparse)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
outdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

load(file.path(outdir, pasteu(exptstr, method, "data.Rda")))
if (method == "scone") {
   mat <- get_normalized(scone_obj3, opt$normalization)
   rm(scone_obj3)
   gc()
}
seed <- 2782472
cl <- clusterMany(mat, isCount = TRUE, reduceMethod = "PCA", nReducedDims = 30,
                  alphas = 0.1,
                  sequential = TRUE, subsample = TRUE, minSizes = 5,
                  clusterFun = c("hierarchical01"), subsampleArgs = list(largeDataset = TRUE),
                  ks = 4:12, random.seed = seed,
                  ncores = ncores)

save(cl, file = file.path(outdir, pasteu(exptstr, "_scone_rsec.Rda")))