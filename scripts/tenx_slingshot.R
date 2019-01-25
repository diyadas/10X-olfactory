# Slingshot lineages on 10X data
# Authors: Diya Das
# Last revised: Wed Jan 23 12:20:31 2019

# Load command-line arguments
library(Seurat)
library(clusterExperiment)
library(slingshot)
library(BiocParallel)
library(optparse)
library(proftools)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--norm", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--clusmethod", type = "character", default = "snn",
              help = "clustering method - snn or rsec")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

datfile <- last_datfile(datdir, exptstr, opt$method, opt$norm, opt$clusmethod)
load(datfile)

reducedDim <- seu@dr$pca@cell.embeddings
clusterLabels <- seu@meta.data$res.0.5
names(clusterLabels) <- rownames(reducedDim)

for (k in 3:20){
    print(paste("Running slingshot for", k, "PCs:"))
    pd1 <- profileExpr(slingOut <- getLineages(reducedDim[,1:k], clusterLabels = clusterLabels,
                      start.clus = 11, end.clus = c(9, 3)))
    print(hotPaths(pd1), value = 'time')
    print(slingLineages(slingOut))
}

