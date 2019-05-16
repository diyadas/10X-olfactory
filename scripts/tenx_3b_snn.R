# Clustering of 10X data using SNN (Seurat)
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Oct 15 18:16:05 2018

# Load command-line arguments
library(scone)
library(zinbwave)
library(Seurat) #packageVersion("Seurat") on bridges is 3.0
library(BiocParallel)
library(optparse)
library(mgcv)

option_list <- list(
  make_option("--expt", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--idfilt", default = FALSE, type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores

print(opt)

if (opt$idfilt) {
  idfiltstr <- "idfiltyes"
} else {
  idfiltstr <- "idfiltno"
}

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

# exclude late-traced cells, different experiment
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "1_se_filtqc", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

#ob_1_se_filtqc_idfiltno.Rda
samples <- colnames(se_filtered)[grep("late", colData(se_filtered)$expt, invert = TRUE)] 
print(head(samples))
print(length(samples))

if (method == "scone") {
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "scone", opt$normalization, "data", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
} else if (method == "zinb"){
#  load(file.path(datdir, pasteu(exptstr, method, "data.Rda")))
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, idfiltstr, "zinb",  "data"),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
}

#load(file.path(datdir, pasteu0(exptstr, "scone", opt$normalization, "data", idfiltstr, ".Rda")))
mat <- get_normalized(scone_obj, opt$normalization, log = FALSE)
mat <- log2(mat + 1)
#rownames(mat)  <- unique(rownames(mat))
rownames(mat) <- make.names(rownames(mat), unique= TRUE)
#mat <- uniquecombs(mat,ordered=FALSE)
rm(scone_obj)
print(dim(mat))
print(head(colnames(mat)))
print(length(setdiff(samples, colnames(mat))))

mat <- mat[, samples]

seed <- 2782472
resolution <- seq(0, 2, 0.1)
#resolution <- (2)

#seu <- CreateSeuratObject(counts= mat, project = exptstr)
seu <- CreateSeuratObject(counts= mat[which(!duplicated(row.names(mat))),], project = exptstr)
#raw.data = mat, min.cells = 1, min.genes = 1,  project = exptstr)
if (method == "scone"){
seu <- FindVariableFeatures(object = seu)
seu <- ScaleData(object = seu, do.center = TRUE, do.scale = FALSE)
seu <- RunPCA(object = seu, seed.use = seed, pc.genes = rownames(seu@data), pcs.compute = 50)
seu <- FindNeighbors(object = seu)
seu <- FindClusters(object = seu, reduction = "pca", 
                    dims.use = 1:50, #this should match K
                    resolution = resolution, print.output = 0, save.SNN = TRUE)
} else if (method == "zinb") {
  seu <- ScaleData(object = seu, do.center = TRUE, do.scale = FALSE)
  seu <- SetDimReduction(object = seu, reduction.type = "zinbwave", 
                         slot = "cell.embeddings",
                         new.data = reducedDim(zinb_obj, "zinbwave")[samples,])
  seu <- SetDimReduction(object = seu, reduction.type = "zinbwave", slot = "key",
                         new.data = "zinbwave")
  seu <- FindClusters(object = seu, reduction.type = "zinbwave", 
                      dims.use = 1:20, #this should match K
                      resolution = resolution, print.output = 0, save.SNN = TRUE)
}


save(seu, file = file.path(datdir, 
	         	   pasteu0(exptstr, opt$method, opt$normalization, "snn",
			   idfiltstr, 
			   format(Sys.time(), "%Y%m%d_%H%M%S"), ".Rda")))