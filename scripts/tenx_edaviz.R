# Visualization of 10X data, pre-ID filtering
# Authors: Diya Das and Rebecca Chance
# Last revised: Mon Apr  8 15:07:07 2019

# Load command-line arguments
library(scone)
library(BiocParallel)
library(optparse)
library(Rtsne)
library(clusterExperiment)

option_list <- list(
  make_option("--expt", default = "", type = "character", help = "Experiment ID"),
  make_option("--normalization", type = "character", help = "name of normalization"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--markerfile", type = "character", help = "marker gene list"),
  make_option("--idfilt", default = FALSE, type = "logical", 
              help = "logical, has sample ID filtering been performed?")
  )

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
ncores <- opt$ncores

if (opt$idfilt) {
  idfiltstr <- ""
} else {
  idfiltstr <- "idfiltno"
}

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

datfiles <<- list.files(path = datdir, 
                        pattern = pasteu(exptstr, "scone", opt$normalization), 
                        "data", idfiltstr,
                        full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
print(opt)
print(dim(scone_obj))

seed <- 2782472

# Plot marker gene heatmap
massivePalette <- massivePalette[-3]

mycoldata <- data.frame(colData(scone_obj), samples = colnames(scone_obj)) %>% arrange(expt, batch)
scone_obj <- scone_obj[, mycoldata$scone_obj]

markers <- intersect(unlist(read.table(file.path("../ref", opt$markerfile))),
                     rownames(scone_obj))

logcounts <- log2(scone_obj@assays$data[[1]] + 1)
breakv <- c(min(logcounts),
            seq(0, quantile(logcounts[logcounts > 0], .99, na.rm = TRUE), length = 50),
            max(logcounts))
breakv <- unique(breakv)

pdf(file = file.path(vizdir,
                     pasteu0(exptstr, "EDA", "markerhm", 
                             "scone", opt$normalization,
                             format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")),
    width = 8.5, height = 11)
plotHeatmap(logcounts,
            clusterFeaturesData = markers, clusterFeatures = FALSE,
            clusterSamples = FALSE, breaks = breakv,
            annLegend = TRUE, overRideClusterLimit = TRUE,
            clusterLegend = list(expt = massivePalette[1:nlevels(mycoldata$expt)], 
                                 batch = massivePalette[1:nlevels(mycoldata$batch)]),
            colData = mycoldata[, c("expt", "batch")],
            labCol = rep("", ncol(scone_obj)))
dev.off()

# t-SNE colored by cluster and time point
rtsne_fx <- function(logcounts, ngenes, perp) {
  set.seed(9887)
  genes.use <- names(vars)[1:ngenes]
  var_data <- logcounts[genes.use, ]
  tsne_data <- Rtsne(t(var_data), 
                       perplexity = perp, max_iter = 10000, num_threads = 0)
  return(tsne_data)
}

vars <- apply(logcounts, 1, var)
vars <- sort(vars, decreasing = TRUE)

tsne_data <- rtsne_fx(logcounts, ngenes = 500, perp = 80)

pdf(file = file.path(vizdir, pasteu0(exptstr, "EDA", "tsne", ngenes, perp,
                                       "scone", opt$normalization, 
                                       format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                       ".pdf")))
if (nlevels(expt) > 1) {
    plot(tsne_data$Y, pch = 19, cex = 0.4, 
         col = alpha(massivePalette[mycoldata$expt], 0.3), 
         xlab = "TSNE 1", ylab = "TSNE 2", 
         main = paste("expt:", ngenes, "genes, perplexity =", perp))
    legend("bottomleft", legend = levels(expt), fill = massivePalette, cex = 0.6)
  }
  
  if (nlevels(batch) > 1) {
    plot(tsne_data$Y, pch = 19, cex = 0.4, 
         col = alpha(massivePalette[mycoldata$batch], 0.3), 
         xlab = "TSNE 1", ylab = "TSNE 2", 
         main = paste("batch:", ngenes, "genes, perplexity =", perp))
    legend("bottomleft", legend = levels(batch), fill = massivePalette, cex = 0.6)
  }
  
  tsne_expr <- data.frame(tsne_data$Y, t(logcounts[markers, ]))
  par(mar = c(2, 2, 1, 1), mfrow = c(1, 1))
  for (gene in markers) {
    p <- ggplot(tsne_expr, aes_string("X1", "X2", colour = gene)) + geom_point(cex= 0.5)
    print(p + t1 +
            scale_colour_gradient2(low = "#053061", mid = "grey95", high = "#67001F") +
            ggtitle(gene))
  }
  dev.off()

