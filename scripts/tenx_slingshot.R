# Differential gene expression of 10X data
# Authors: Diya Das
# Last revised: Mon Oct 15 22:19:11 2018

# Load command-line arguments
library(Seurat)
library(clusterExperiment)
library(slingshot)
library(BiocParallel)
library(optparse)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double")
)

print(sessionInfo())

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
outdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

load(file.path(outdir, pasteu(exptstr, opt$method, "snn.Rda")))
sce <- get(load(file.path(outdir, pasteu(exptstr, opt$method, "sce_res05.Rda"))))
reducedDim(sce) <- seu@dr$pca@cell.embeddings[,1:6]
slingOut <- slingshot(sce, "res.0.5", start.clus = 12, end.clus = c(9,2))
save(slingOut, file = file.path(outdir, 
                                    pasteu(exptstr, method, 
                                           opt$normalization,"res05","slingshot", format(Sys.time(), "%Y%m%d_%H%M%S"),".Rda")))


