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
load(file.path(datdir, pasteu(exptstr, opt$method, "snn.Rda")))
load(file.path(datdir, pasteu(exptstr, "se_filtered.Rda")))


names(batch) <- colnames(seu@data)
names(expt) <- colnames(seu@data)
clus.labels <- sort(clus.labels)
batch <- batch[names(clus.labels)]
expt <- expt[names(clus.labels)]

seed <- 2782472

# heatmaps

oe_markers <- intersect(unlist(read.table("../ref/oe_markers32+regen.txt")), 
                        rownames(seu@data))
breakv <- c(min(seu@data), 
            seq(0, quantile(seu@data[seu@data > 0], .99, na.rm = TRUE), length = 50), 
            max(seu@data))

pdf(file = file.path(vizdir, pasteu(exptstr, "marker", Sys.Date(), ".pdf")), 
    width = 6, height = 5.5)
plotHeatmap(seu@data[oe_markers, names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels, expt = expt, batch = batch), clusterLegend = list(cluster = bigPalette, expt = cole), annLegend = TRUE)
dev.off()

# t-SNE colored by cluster and time point
rtsne_fx <- function(cmobj, clusters, ngenes = 500) {
  set.seed(9887)
  vars <- apply(seu@data[,names(clusters)], 1, var)
  vars <- sort(vars, decreasing = TRUE)
  var_data <- seu@data[names(vars)[1:ngenes],]
  tsne_data <- Rtsne(t(var_data[,names(clusters)]), 
                     perplexity = 10, max_iter = 1000)
  return(tsne_data)
}

tsne_data <- rtsne_fx(seu, clus.labels, ngenes = 500)
pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "clus", method, opt$norm, Sys.Date(), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = bigPalette[clus.labels], xlab = "TSNE 1", ylab = "TSNE 2")
legend("topleft", legend = levels(clus.labels), fill = bigPalette, cex = 0.5)
dev.off()

pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "batch", method, opt$norm, Sys.Date(), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = colb[batch], xlab = "TSNE 1", ylab =" TSNE 2")
legend("bottomleft", legend = levels(batch), fill = colb, cex = 0.6)
dev.off()

pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "expt", method, opt$norm, Sys.Date(), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = cole[expt], xlab = "TSNE 1", ylab =" TSNE 2")
legend("bottomleft", legend = levels(expt), fill = cole, cex = 0.6)
dev.off()
