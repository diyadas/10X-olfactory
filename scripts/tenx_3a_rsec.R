# Clustering of 10X data using RSEC
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Jul  9 16:07:04 2018

# Load command-line arguments
library(scone)
library(SingleCellExperiment)
library(clusterExperiment)
library(BiocParallel)
library(optparse)

print(sessionInfo())

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--sequential", default = FALSE, type = "logical", help = "TRUE or FALSE"),
  make_option("--subsample", default = FALSE, type = "logical", help = "TRUE or FALSE"),
  make_option("--idfilt", type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

if (opt$idfilt) {
  idfiltstr <- "idfiltyes"
} else {
  idfiltstr <- "idfiltno"
}

print(opt)
mytimestamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz="America/Los_Angeles")
print(paste("Files produced by this script will be timestamped:", mytimestamp))

# exclude late-traced cells, different experiment

datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "1_se_filtqc", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
samples <- colnames(se_filtered)
samples <- colnames(se_filtered)[grep("late", colData(se_filtered)$expt, invert = TRUE)] 
se_filtered <- se_filtered[, samples]

datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, "scone", opt$normalization, "data", idfiltstr),
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

mat <- get_normalized(scone_obj, opt$normalization, log = FALSE)
rm(scone_obj)
gc()
mat <- mat[, samples]
rownames(mat) <- make.names(rownames(mat), unique= TRUE)
reduceMethod = "PCA"

sce <- SingleCellExperiment(assays = list(counts = mat))

if (method == "zinb") {

  datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, idfiltstr, method, "data"),
                        full.names = TRUE)
  datfile <- datfiles[length(datfiles)]
  print(paste("Loading this data file: ", datfile))
  load(datfile)

  reducedDims(sce) <- SimpleList(zinb = reducedDim(zinb_obj, "zinbwave")[samples,])
  reduceMethod = "zinb"
}

colData(sce) <- DataFrame(expt = colData(se_filtered)$expt,
                          batch = colData(se_filtered)$batch,
                          samples = samples)

cl <- RSEC(sce, k0s = seq(10, 25, by = 3), alphas = c(0.1, 0.2, 0.3),
           reduceMethod = reduceMethod,
           transFun = function(x) log2(x + 1), # should override isCount=FALSE, or could write isCount=TRUE
           nReducedDims = 20,
#           nReducedDims = 30,
           sequential = opt$sequential, subsample = opt$subsample,
           subsampleArgs = list(resamp.num = 50, clusterFunction = "kmeans"),
           betas = c(0.8),
           clusterFunction = c("hierarchical01", "pam"),
           minSizes = 10,
           ncores = ncores,
           consensusProportion = 0.7, consensusMinSize = 10, # with sequential turned off, then consensus might result in many -1s
           dendroReduce = "PCA", dendroNDims = 30,
	   mergeMethod = "none",
       #    mergeMethod = "adjP", 
       #    mergeCutoff = 0.1,
       #    mergeLogFCcutoff = 1, 
       #    mergeDEMethod = "limma-voom",
	   random.seed = 2357891,
           verbose = TRUE, run = TRUE)
cl
print(colnames(cl@clusterMatrix))
subsamplestr <- ifelse (opt$subsample, "sub", "nosub")
seqstr <- ifelse(opt$sequential, "seq", "noseq")

traceback()

save(cl, subsamplestr, seqstr, 
     file = file.path(datdir, pasteu0(exptstr, opt$method, opt$normalization,
                                      "rsec",
                                      idfiltstr,
				      mytimestamp,
                                      ".Rda")))
