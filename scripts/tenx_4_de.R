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
  idfiltstr <- ""
} else {
  idfiltstr <- "idfiltno"
}

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")
datfiles <- list.files(path = datdir, pattern = pasteu(exptstr, method, opt$normalization, "snn", idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

load(file.path(datdir, pasteu(exptstr, "se_filtered.Rda")))
se_filtered <- se_filtered[, colnames(seu@data)]
metadata <- data.frame(seu@meta.data, expt = colData(se_filtered)$expt,
                                      batch = colData(se_filtered)$batch,
                      samples = colnames(seu@data), row.names = "samples")
metadata <- metadata[colnames(seu@raw.data),]

#ce <- ClusterExperiment(seu@raw.data, clusters = seu@meta.data[,"res.0.5"], 
#                        transformation = function(x) log2(x + 1))
ce <- ClusterExperiment(seu@raw.data, clusters = seu@meta.data[,"res.2"],
                        transformation = function(x) log2(x + 1))
colData(ce) <- DataFrame(metadata)

ce <- makeDendrogram(ce, reduceMethod = "var", nDims = 1000)

de_ce <- getBestFeatures(ce, contrastType = "OneAgainstAll", whichAssay = 1, 
                         DEMethod = "limma", number = 100)

rownames(se_filtered) <- rowData(se_filtered)$Symbol
colData(ce)$batch <- colData(se_filtered)$batch

save(ce, file = file.path(datdir, pasteu0(exptstr, method,
                                           opt$normalization,"res05","cmobj", format(Sys.time(), "%Y%m%d_%H%M%S"),".Rda")))

write.table(de_ce, file = file.path(datdir, 
                                    pasteu0(exptstr, method, 
                                           opt$normalization,"res05","DE_OneAgainstAll", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

#de_ce <- getBestFeatures(ce, contrastType = "Pairs", whichAssay = 1, 
#                         DEMethod = "limma", number = 100)
#write.table(de_ce, file = file.path(datdir, 
#                                    pasteu(exptstr, method, 
#                                           opt$normalization,"DE_Pairs", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
#            sep = "\t", quote = FALSE, row.names = FALSE)

