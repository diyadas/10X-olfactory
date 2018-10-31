# Visualization of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Wed Oct 31 09:35:01 2018

# Load command-line arguments
library(scone)
library(clusterExperiment)
library(BiocParallel)
library(optparse)
library(Rtsne)
library(Seurat)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--norm", type = "character", help = "name of normalization"),
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

metadata <- data.frame(seu@meta.data, expt = expt, batch = batch,  
	    	      samples = colnames(seu@data), row.names = "samples")
metadata <- metadata[with(metadata, order(res.0.5, expt, batch)),]
metadata <- metadata[, c('res.0.5', 'expt', 'batch')]

seed <- 2782472

# heatmaps

oe_markers <- intersect(unlist(read.table("../ref/oe_markers32+regen.txt")), 
                        rownames(seu@data))
breakv <- c(min(seu@data), 
            seq(0, quantile(seu@data[seu@data > 0], .99, na.rm = TRUE), length = 50), 
            max(seu@data))
breakv <- unique(breakv)
pdf(file = file.path(vizdir, pasteu(exptstr, "markerhm", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")), 
    width = 6, height = 5.5)
plotHeatmap(seu@data[oe_markers, rownames(metadata)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = metadata, 
				 clusterLegend = list(res.0.5 = bigPalette, expt = cole, batch = bigPalette), annLegend = TRUE)
dev.off()

# t-SNE colored by cluster and time point


if (is.null(seu@dr$tsne)){
ngenes = 500
vars <- apply(seu@data, 1, var)
vars <- sort(vars, decreasing = TRUE)
genes.use <- names(vars)[1:ngenes]

seu <- RunTSNE(seu, reduction.use = "pca", dims.use = 1:50,
  genes.use = genes.use, seed.use = seed, tsne.method = "Rtsne", perplexity = 10, max_iter = 10000,
  dim.embed = 2, reduction.name = "tsne")

save(seu, file.path(datdir, pasteu(exptstr, opt$method, "snn.Rda")))
}

pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "clus", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = bigPalette[factor(seu@meta.data[,"res.0.5"])], xlab = "TSNE 1", ylab = "TSNE 2")
legend("topleft", legend = levels(factor(seu@meta.data[,"res.0.5"])), fill = bigPalette, cex = 0.5)
dev.off()

pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "batch", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = colb[batch], xlab = "TSNE 1", ylab =" TSNE 2")
legend("bottomleft", legend = levels(batch), fill = colb, cex = 0.6)
dev.off()

pdf(file = file.path(vizdir, pasteu(exptstr, "tsne", "expt", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
plot(tsne_data$Y, pch = 19, cex = 0.4, col = cole[expt], xlab = "TSNE 1", ylab =" TSNE 2")
legend("bottomleft", legend = levels(expt), fill = cole, cex = 0.6)
dev.off()
