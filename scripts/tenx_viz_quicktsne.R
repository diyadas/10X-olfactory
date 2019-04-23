# Visualization of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Thu Feb 28 15:58:37 2019

# Load command-line arguments
library(scone)
library(BiocParallel)
library(optparse)
library(Rtsne)
library(Seurat)
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

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

datfiles <<- list.files(path = datdir, pattern = glob2rx(pasteu(exptstr, method, opt$norm, 
                                                        opt$clusmethod, "*viz*")), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

print(opt)

seed <- 2782472

# Plot marker gene heatmap
massivePalette <- massivePalette[-3]

## Set cluster legend colors
clusterLegend <- setNames(rep(list(massivePalette), ncol(colData(cl))), colnames(colData(cl)))
#if (exptstr == "ob" & length(unique(colData(cl)[, 1])) < length(colRKC)) {
#  clusterLegend <- setNames(rep(list(colRKC), ncol(colData(cl))), colnames(colData(cl)))
#}

# clusterLegend[["expt"]] <- cole

markers <- intersect(unlist(read.table(file.path("../ref", opt$markerfile))),
                     rownames(cl))
dat <- transformData(cl)
breakv <- c(min(dat),
            seq(0, quantile(dat[dat > 0], .99, na.rm = TRUE), length = 50),
            max(dat))
breakv <- unique(breakv)

print(length(massivePalette))
print(length(unique(primaryCluster(cl))))
pdf(file = file.path(vizdir,
                     pasteu0(exptstr, "markerhm", method, opt$norm,
                             opt$clusmethod,
                             format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")),
    width = 8.5, height = 11)
plotHeatmap(cl,
            whichClusters = "all",
            clusterSamplesData = opt$samplesort,
            clusterFeaturesData = markers, clusterFeatures = FALSE,
            breaks = breakv, 
            clusterLegend = clusterLegend,
            annLegend = TRUE, overRideClusterLimit = TRUE,
            colData = which(colnames(colData(cl)) %in% c("expt", "batch")),
            labCol = rep("", ncol(cl)))
dev.off()

# t-SNE colored by cluster and time point

#ngenesvec <- c(500, 1000, 5000, nrow(cl))
#perpvec <- seq(50, 100, 25)

ngenesvec = 500
perpvec  = 75

vars <- apply(transformData(cl), 1, var)
vars <- sort(vars, decreasing = TRUE)

numclus <- length(setdiff(unique(primaryCluster(cl)), "-1"))

params <- expand.grid(perp = perpvec, ngenes = ngenesvec)
lapply(1:nrow(params), function(x) {
  message(params[x,])
  ngenes <- params[x, "ngenes"]
  perp <- params[x, "perp"]

  pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", ngenes, perp,
                                       method, opt$clusmethod, opt$norm, 
                                       format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                       ".pdf")))
  
  plot(cl@reducedDims$tsne[[x]]$Y, pch = 19, cex = 0.4, 
       col = alpha(massivePalette[factor(primaryCluster(cl))], 0.3), 
       xlab = "TSNE 1", ylab = "TSNE 2", 
       main = paste(numclus, "clusters:", ngenes, "genes, perplexity =", perp))
  plot.new()
  showPalette(massivePalette[1:100],cex=0.6)
  showPalette(massivePalette[101:200],cex=0.6)
  showPalette(massivePalette[201:300],cex=0.6)
  showPalette(massivePalette[301:400],cex=0.6)
  showPalette(massivePalette[401:length(massivePalette)],cex=0.6)
  plot.new()
  legend("topleft", legend = levels(factor(primaryCluster(cl))), 
         fill = massivePalette, cex = 0.5, ncol = 5)
  batch <- factor(colData(cl)$batch)
  expt <- factor(colData(cl)$expt)
  
  if (nlevels(expt) > 1) {
    plot(cl@reducedDims$tsne[[x]]$Y, pch = 19, cex = 0.4, 
        col = alpha(massivePalette[expt], 0.3), 
        xlab = "TSNE 1", ylab = "TSNE 2", 
        main = paste("expt:", ngenes, "genes, perplexity =", perp))
   legend("bottomleft", legend = levels(expt), fill = massivePalette, cex = 0.6)
 }
 
 if (nlevels(batch) > 1) {
   plot(cl@reducedDims$tsne[[x]]$Y, pch = 19, cex = 0.4, 
        col = alpha(massivePalette[batch], 0.3), 
        xlab = "TSNE 1", ylab = "TSNE 2", 
        main = paste("batch:", ngenes, "genes, perplexity =", perp))
   legend("bottomleft", legend = levels(batch), fill = massivePalette, cex = 0.6)
 }
 
 dat <- data.frame(cl@reducedDims$tsne[[x]]$Y, t(transformData(cl)[markers, ]))
 par(mar = c(2, 2, 1, 1), mfrow = c(1, 1))
 for (gene in markers) {
   p <- ggplot(dat, aes_string("X1", "X2", colour = gene)) + geom_point(cex= 0.5)
   print(p + t1 +
           scale_colour_gradient2(low = "#053061", mid = "grey95", high = "#67001F") +
           ggtitle(gene))
 }
  dev.off()
  
})

