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
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")
load(file.path(outdir, pasteu(exptstr, opt$method, "snn.Rda")))
load(file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
#load(file.path(outdir, pasteu(exptstr, method, "data.Rda")))

clus.labels <- sort(clus.labels)
batch <- batch[clus.labels]
expt <- expt[clus.labels]

seed <- 2782472

# heatmaps

oe_markers <- intersect(unlist(read.table("../ref/oe_markers32+regen.txt")),rownames(seu@data))
breakv <- c(min(seu@data), seq(0, quantile(seu@data[seu@data > 0], .99, na.rm = TRUE), length = 50), max(seu@data))

pdf(file = file.path(vizdir, pasteu(exptstr, "selectmarker", Sys.Date(), ".pdf")), width = 6, height = 5.5)
plotHeatmap(seu@data[oe_markers, names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels, expt = expt, batch = batch), clusterLegend = list(expt = cole), annLegend = TRUE)
dev.off()

# t-SNE colored by cluster and time point
