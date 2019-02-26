# Clustering of 10X data using RSEC
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Mon Jul  9 16:07:04 2018

# Load command-line arguments
library(scone)
library(clusterExperiment)
library(BiocParallel)
library(optparse)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--sequential", default = FALSE, type = "logical", help = "TRUE or FALSE"),
  make_option("--subsample", default = FALSE, type = "logical", help = "TRUE or FALSE")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
outdir <- file.path("../output", exptstr, "data")
ncores <- opt$ncores

#register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

message(str(args))

load(file.path(outdir, pasteu(exptstr, method, opt$normalization, "data.Rda")))
if (method == "scone") {
   mat <- get_normalized(scone_obj, opt$normalization)
   rm(scone_obj)
   gc()
}
seed <- 2782472

cl <- RSEC(mat, k0s = seq(10, 25, by=3), alphas = c(0.1, 0.2, 0.3),
                reduceMethod = "PCA",
                transFun = function(x) log2(x+1), # should override isCount=FALSE, or could write isCount=TRUE
                nReducedDims = 30,
                sequential = opt$sequential, subsample = opt$subsample,
                subsampleArgs = list(resamp.num = 50, clusterFunction="kmeans"),
                betas = c(0.8),
                clusterFunction = c("hierarchical01", "tight"),  
		minSizes = 1,
                ncores = ncores,
                consensusProportion = 0.7, consensusMinSize = 10, # with sequential turned off, then consensus might result in many -1s
                dendroReduce = "PCA", dendroNDims = 30,
                mergeMethod = "adjP", mergeCutoff = 0.1,
                mergeLogFCcutoff = 1, random.seed = 2357891,
                mergeDEMethod = "limma-voom", 
                verbose = TRUE, run = TRUE)

subsamplestr <- ifelse (opt$subsample, "sub", "nosub")
seqstr <- ifelse(opt$sequential, "seq", "noseq")

load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
se_filtered <- se_filtered[, colnames(cl)]
colData(cl) <- DataFrame(cl@clusterMatrix, colData(cl),
			 expt = colData(se_filtered)$expt,
                         batch = colData(se_filtered)$batch,
                         samples = colnames(cl))

save(cl, subsamplestr, seqstr, file = file.path(outdir, pasteu0(exptstr, opt$method, opt$normalization, "rsec",
                                               format(Sys.time(), "%Y%m%d_%H%M%S") ,".Rda")))