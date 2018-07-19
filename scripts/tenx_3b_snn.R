# Clustering of 10X data using SNN (Seurat)
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Thu Jul 19 13:46:01 2018

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
  rm(scone_obj)
}
seed <- 2782472

seu <- CreateSeuratObject(raw.data = mat, min.cells = 1, min.genes = 1, 
                          project = exptstr)

seu <- ScaleData(seu)
seu <- RunPCA(seu, seed.use = seed, pc.genes = rownames(seu@data))
resolution = 1
seu <- FindClusters(object = seu, reduction.type = "pca", 
                    dims.use = 1:20, #this should match K
                    resolution = 1, print.output = 0, save.SNN = TRUE)
clus.labels <- seu@ident
names(clus.labels) <- paste0(names(seu@ident))

save(seu, clus.labels, file = file.path(outdir, pasteu(exptstr, opt$method, "snn.Rda")))