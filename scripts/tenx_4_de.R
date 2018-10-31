# Differential gene expression of 10X data
# Authors: Diya Das
# Last revised: Mon Oct 15 22:19:11 2018

# Load command-line arguments
library(Seurat)
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

load(file.path(outdir, pasteu(exptstr, opt$method, "snn.Rda")))

ce <- ClusterExperiment(seu@raw.data, clusters = seu@meta.data[,"res.0.5"], 
                        transformation = function(x) log2(x + 1))
ce <- makeDendrogram(ce, reduceMethod = "var", nDims = 1000)

de_ce <- getBestFeatures(ce, contrastType = "OneAgainstAll", whichAssay = 1, 
                         DEMethod = "limma", number = 100)
write.table(de_ce, file = file.path(outdir, 
                                    pasteu(exptstr, method, 
                                           opt$normalization,"res05","DE_OneAgainstAll", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)


#de_ce <- getBestFeatures(ce, contrastType = "Pairs", whichAssay = 1, 
#                         DEMethod = "limma", number = 100)
#write.table(de_ce, file = file.path(outdir, 
#                                    pasteu(exptstr, method, 
#                                           opt$normalization,"DE_Pairs", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
#            sep = "\t", quote = FALSE, row.names = FALSE)

