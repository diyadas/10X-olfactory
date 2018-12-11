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

print(paste("This job is using", ncores, "cores of a node."))

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

# exclude late-traced cells, different experiment
load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
samples <- colnames(se_filtered)[grep("late", colData(se_filtered)$expt, invert = TRUE)] 

load(file.path(outdir, pasteu(exptstr, "scone", "data.Rda")))
mat <- get_normalized(scone_obj3, opt$normalization)
mat <- log2(mat + 1)
rm(scone_obj3)

if (method == "zinb"){
  load(file.path(outdir, pasteu(exptstr, method, "data.Rda")))
}

mat <- mat[, samples]

seed <- 2782472
resolution <- seq(0, 2, 0.1)

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
                         new.data = reducedDim(zinb_obj, "zinbwave")[samples,])
  seu <- SetDimReduction(object = seu, reduction.type = "zinbwave", slot = "key",
                         new.data = "zinbwave")
  seu <- FindClusters(object = seu, reduction.type = "zinbwave", 
                      dims.use = 1:20, #this should match K
                      resolution = resolution, print.output = 0, save.SNN = TRUE)
}

clus.labels <- seu@ident
names(clus.labels) <- paste0(names(seu@ident))

save(seu, clus.labels, file = file.path(outdir, 
	  	       	      		pasteu(exptstr, opt$method, opt$normalization, "snn", 
					       format(Sys.time(), "%Y%m%d_%H%M%S") ,".Rda")))