# Differential gene expression of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Thu May 16 2019

# Load command-line arguments
library(Seurat)
library(clusterExperiment)
library(BiocParallel)
library(optparse)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--clusmethod", type = "character", default = "snn",
              help = "clustering method - snn or rsec"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--whichmerge", type = "character", default = "adjP_cutoff_0.01", help = "choice post rsec"),
  make_option("--idfilt", default = FALSE, type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores
clusmethod <- opt$clusmethod
whichmerge <- opt$whichmerge

print(opt)

if (opt$idfilt) {
  idfiltstr <- "idfiltyes"
} else {
  idfiltstr <- "idfiltno"
}

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")
if (opt$normalization == "scone") {
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, method, opt$normalization, 
                                                        opt$clusmethod, idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
} else if (opt$normalization == "zinb") {
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, method, "data", idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
}
#load(file.path(datdir, pasteu0(exptstr, "1_se_filtqc", idfiltstr, ".Rda")))

if (opt$clusmethod == "snn") {
datfiles <- list.files(path = datdir, pattern = pasteu(exptstr, method, opt$normalization, "snn", idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
} else if (opt$clusmethod == "rsec") {
   if (opt$method == "scone") {
          #ob_scone_none,fq,qc_k=2,no_bio,batch_rsec_locfdr_mergecutoff_0.5_20190425_155428.Rda #0.08 here
   datfiles <- list.files(path = datdir, 
                       pattern = pasteu(exptstr, method, opt$normalization, clusmethod, whichmerge), full.names = TRUE)
    datfile <- datfiles[length(datfiles)]
    print(paste("Loading this data file: ", datfile))
    load(datfile)
   } else if (opt$method == "zinb") {
        datfiles <- list.files(path = datdir, 
                       pattern = pasteu(exptstr, method, opt$normalization, clusmethod, whichmerge), full.names = TRUE)
        datfile <- datfiles[length(datfiles)]
        print(paste("Loading this data file: ", datfile))
        load(datfile)
   }
}

#pare down to filtered
datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "1_se_filtqc", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

seu@data <- GetAssayData(object = seu)
se_filtered <- se_filtered[, colnames(seu@data)]
seu@meta.data <- GetAssayData(object = seu, slot = "counts")
metadata <- data.frame(seu@meta.data, expt = colData(se_filtered)$expt,
                                      batch = colData(se_filtered)$batch,
                      samples = colnames(seu@data), row.names = "samples")
metadata <- metadata[colnames(seu@raw.data),]

if (opt$clusmethod == "rsec") {
ce <- ClusterExperiment(cl2, clusters = primaryCluster)
ce <- makeDendrogram(ce, reduceMethod = "var", nDims = 1000)
de_onevall <- getBestFeatures(ce, contrastType = "OneAgainstAll", whichAssay = 1,
                         DEMethod = "limma", number = 100)
} else if (opt$clusmethod == "snn") {
#ce <- ClusterExperiment(seu@raw.data, clusters = seu@meta.data[,"res.0.5"], 
#                        transformation = function(x) log2(x + 1))
ce <- ClusterExperiment(seu@raw.data, clusters = seu@meta.data[,"res.2"],
                        transformation = function(x) log2(x + 1))
colData(ce) <- DataFrame(metadata)

ce <- makeDendrogram(ce, reduceMethod = "var", nDims = 1000)

de_ce <- getBestFeatures(ce, contrastType = "OneAgainstAll", whichAssay = 1,
                         DEMethod = "limma", number = 100)
}
rownames(se_filtered) <- rowData(se_filtered)$Symbol
colData(ce)$batch <- colData(se_filtered)$batch


save(ce, file = file.path(datdir, pasteu0(exptstr, method,
                                           opt$normalization,"res05","cmobj", whichmerge,  format(Sys.time(), "%Y%m%d_%H%M%S"),".Rda")))

write.table(de_ce, file = file.path(datdir, 
                                    pasteu0(exptstr, method, 
                                           opt$normalization,"res05","DE_OneAgainstAll", whichmerge,  format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

#de_ce <- getBestFeatures(ce, contrastType = "Pairs", whichAssay = 1, 
#                         DEMethod = "limma", number = 100)
#write.table(de_ce, file = file.path(datdir, 
#                                    pasteu(exptstr, method, 
#                                           opt$normalization,"DE_Pairs", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt")),
#            sep = "\t", quote = FALSE, row.names = FALSE)

