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
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--method", type = "character", help = "scone or zinb"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--markerfile", default = "oe_markers32+regen.txt", type = "character",
              help = "marker gene list"),
  make_option("--clusmethod", type = "character", default = "snn",
              help = "clustering method - snn or rsec"),
  make_option("--seures", type = "character", default = "res.0.5",
              help = "Seurat, which resolution to use for primary clustering"),
  make_option("--samplesort", type = "character", 
              help = "argument clusterSamplesData, i.e. dendrogramValue or primaryCluster"),
  make_option("--idfilt", default = FALSE, type = "logical", help = "logical, has sample ID filtering been performed?")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (opt$idfilt) {
  idfiltstr <- "idfiltyes"
} else {
  idfiltstr <- "idfiltno"
}

exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

#chosen normalization file
datfiles <<- list.files(path = datdir, 
                        pattern = pasteu(exptstr, opt$method, opt$normalization, 
                                         opt$clusmethod, 
                                         idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)

#
datfiles <<- list.files(path = datdir,
                        pattern = pasteu(exptstr, "1_se_filtqc", idfiltstr), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
#load(file.path(datdir, pasteu0(exptstr, "1_se_filtqc", idfiltstr, ".Rda")))

if (opt$clusmethod == "snn"){
  seu@data <- GetAssayData(seu = seu)
  se_filtered <- se_filtered[, colnames(seu@data)]
  metadata <- seu[, c(grep("^res", colnames(seu))), drop = FALSE]
  #metadata <- seu@meta.data[, c(grep("^res", colnames(seu@meta.data))), drop = FALSE]
  if (opt$seures %in% colnames(metadata)) {
     seures <- 0.5
  } else {
    message("WARNING: seures not in range of clusterings, using seures as below")
    seures <- colnames(metadata)[1]
    message(paste("seures:", seures))
  }
  counts <- 2 ^ (seu@data) - 1
  cl <- ClusterExperiment(counts,
                          clusters = as.matrix(metadata),
                          primaryIndex = which(colnames(metadata) == seures),
                          transformation = function(x) log2(x + 1))
  metadata <- data.frame(expt = colData(se_filtered)$expt, 
                         batch = colData(se_filtered)$batch,  
                         samples = colnames(seu@data))
  colData(cl) <- DataFrame(metadata)
if (opt$method == "zinb"){ reducedDim(cl, "zinbwave") <- seu@dr$zinbwave@cell.embeddings
}}


seed <- 2782472

# Plot marker gene heatmap
massivePalette <- massivePalette[-3]
#cl <- cl2
cl <- recolorMassive(cl)
#cl <- recolorMassive_cl(cl)

markers <- intersect(unlist(read.table(file.path("../ref", opt$markerfile))),
                     rownames(cl))
dat <- transformData(cl)
breakv <- c(min(dat),
            seq(0, quantile(dat[dat > 0], .99, na.rm = TRUE), length = 50),
            max(dat))
breakv <- unique(breakv)

pdf(file = file.path(vizdir,
                     pasteu0(exptstr, "markerhm", method, opt$normalization,
                             opt$clusmethod,
                             format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")),
    width = 8.5, height = 11)
plotHeatmap(cl,
            whichClusters = "all",
            clusterSamplesData = opt$samplesort,
            clusterFeaturesData = markers, clusterFeatures = FALSE,
            breaks = breakv, 
            #            clusterLegend = clusterLegend,
            annLegend = TRUE, overRideClusterLimit = TRUE,
            colData = which(colnames(colData(cl)) %in% c("expt", "batch")),
            labCol = rep("", ncol(cl)))
dev.off()

# t-SNE colored by cluster and time point
rtsne_fx <- function(cmobj, ngenes, perp) {
  set.seed(9887)
  genes.use <- names(vars)[1:ngenes]
  if (opt$method == "scone") {
    var_data <- transformData(cmobj)[genes.use, ]
    tsne_data <- Rtsne(t(var_data), 
                       perplexity = perp, max_iter = 10000, num_threads = 0)
  } else if (opt$method == "zinb") {
    W <- cl@reducedDims$zinbwave
    tsne_data <- Rtsne(W, pca = FALSE,
                       perplexity = perp, max_iter = 10000, num_threads = 0)
  }
  return(tsne_data)
}

if (!is.null(cl@reducedDims$tsne)) {
  print("recomputing t-SNE")
}
cl@reducedDims$tsne <- list()

ngenesvec <- c(500, 1000, 5000, nrow(cl))
perpvec <- seq(50, 100, 25)
vars <- apply(transformData(cl), 1, var)
vars <- sort(vars, decreasing = TRUE)

params <- expand.grid(perp = perpvec, ngenes = ngenesvec)
lapply(1:nrow(params), function(x) {
  message(params[x,])
  ngenes <- params[x, "ngenes"]
  perp <- params[x, "perp"]
  cl@reducedDims$tsne[[x]] <- rtsne_fx(cl, ngenes, perp)
  
  pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", ngenes, perp,
                                       method, opt$clusmethod, opt$normalization, 
                                       format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                       ".pdf")))
  
  plot(cl@reducedDims$tsne[[x]]$Y, pch = 19, cex = 0.4, 
       col = alpha(massivePalette[factor(primaryCluster(cl))], 0.3), 
       xlab = "TSNE 1", ylab = "TSNE 2", 
       main = paste("cluster:", ngenes, "genes, perplexity =", perp))
  legend("topleft", legend = levels(factor(primaryCluster(cl))), 
         fill = massivePalette, cex = 0.5)
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

save(cl, file = gsub(".Rda", "_viz.Rda", datfile))