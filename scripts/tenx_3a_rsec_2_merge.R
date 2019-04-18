# Visualization of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Thu Feb 28 15:58:37 2019

# Load command-line arguments
library(scone)
library(BiocParallel)
library(optparse)
library(Rtsne)
library(clusterExperiment)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--norm", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--markerfile", default = "oe_markers32+regen.txt", type = "character",
              help = "marker gene list"),
  make_option("--clusmethod", type = "character", default = "snn",
              help = "clustering method - snn or rsec"),
  make_option("--seures", type = "character", default = "res.0.5",
              help = "Seurat, which resolution to use for primary clustering"),
  make_option("--samplesort", type = "character", 
              help = "argument clusterSamplesData, i.e. dendrogramValue or primaryCluster")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
ncores <- opt$ncores
mytimestamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz="America/Los_Angeles")
print(paste("Files produced by this script will be timestamped:", mytimestamp))


register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

datfiles <<- list.files(path = datdir, pattern = glob2rx(pasteu(exptstr, method, opt$norm, 
                                                        opt$clusmethod, "2019*")), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

print(opt)

seed <- 2782472

# Plot marker gene heatmap
massivePalette <- massivePalette[-3]

## Set cluster legend colors


markers <- intersect(unlist(read.table(file.path("../ref", opt$markerfile))),
                     rownames(cl))
dat <- transformData(cl)
breakv <- c(min(dat),
            seq(0, quantile(dat[dat > 0], .99, na.rm = TRUE), length = 50),
            max(dat))
breakv <- unique(breakv)

mergemethod <- "locfdr"

cutoffvec <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
for (cutoff in cutoffvec) {
pdf(file = file.path(vizdir,
                     pasteu0(exptstr, method, opt$norm,
                             opt$clusmethod, "plotdendro", mergemethod, "cutoff", cutoff,
                             mytimestamp, ".pdf")),
    width = 20, height = 20)

print(cutoff)
print("clusterings pre-merging")
print(colnames(cl@clusterMatrix))

cl2 <- mergeClusters(cl, mergeMethod = mergemethod,
      			DEMethod = "limma-voom",
			plotInfo = "mergeMethod",
			logFCcutoff = 1, cutoff = cutoff,
			leafType = "clusters", plotType = "name", 
			ncores = ncores, random.seed  = 2357891)
print("clusterings post-merging")
print(table(primaryCluster(cl2)))
print(colnames(cl2@clusterMatrix))
flush.console() 
print("-----------------------------")
#plotDendrogram(cl2,leafType="clusters",mergeInfo='mergeMethod', plotType="name", use.edge.lengths = FALSE)

cl2 <- recolorMassive(cl2)


plotHeatmap(cl2,
            whichClusters = "all",
            clusterSamplesData = opt$samplesort,
            clusterFeaturesData = markers, clusterFeatures = FALSE,
            breaks = breakv, 
            annLegend = TRUE, overRideClusterLimit = TRUE,
            colData = which(colnames(colData(cl2)) %in% c("expt", "batch")),
            labCol = rep("", ncol(cl2)))
dev.off()

save(cl2, file = file.path(datdir,
                     pasteu0(exptstr, method, opt$norm,
                             opt$clusmethod, mergemethod, "mergecutoff", cutoff,
                             mytimestamp, ".Rda")))
}