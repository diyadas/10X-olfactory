# Clustering of 10X data using SNN (Seurat)
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Oct 15 18:16:05 2018

# Load command-line arguments
library(scone)
library(zinbwave)
library(Seurat)
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

load(file.path(outdir, pasteu(exptstr, "scone", "data.Rda")))
mat <- get_normalized(scone_obj3, opt$normalization)
mat <- log2(mat + 1)
rm(scone_obj3)

if (method == "zinb"){
  load(file.path(outdir, pasteu(exptstr, method, "data.Rda")))
}

seed <- 2782472
resolution <- seq(0, 10, 0.5)

seu <- CreateSeuratObject(raw.data = mat, min.cells = 1, min.genes = 1, 
                          project = exptstr)
if (method == "scone"){
seu <- ScaleData(seu)
seu <- RunPCA(seu, seed.use = seed, pc.genes = rownames(seu@data))
seu <- FindClusters(object = seu, reduction.type = "pca", 
                    dims.use = 1:20, #this should match K
                    resolution = resolution, print.output = 0, save.SNN = TRUE)
} else if (method == "zinb") {
  seu <- SetDimReduction(object = seu, reduction.type = "zinbwave", 
                         slot = "cell.embeddings",
                         new.data = reducedDim(zinb_obj, "zinbwave"))
  seu <- SetDimReduction(object = seu, reduction.type = "zinbwave", slot = "key",
                         new.data = "zinbwave")
  seu <- FindClusters(object = seu, reduction.type = "zinbwave", 
                      dims.use = 1:20, #this should match K
                      resolution = resolution, print.output = 0, save.SNN = TRUE)
}

clus.labels <- seu@ident
names(clus.labels) <- paste0(names(seu@ident))

save(seu, clus.labels, file = file.path(outdir, pasteu(exptstr, opt$method, "snn.Rda")))