# Visualization of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Wed Oct 31 09:35:01 2018

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
  			      help = "clustering method - snn or rsec")
)

opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
method <- opt$method
datdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz")
ncores <- opt$ncores

register(MulticoreParam(workers = opt$ncores))

source("tenx_helper.R")

datfiles <<- list.files(path = datdir, pattern = pasteu(exptstr, method, opt$norm, opt$clusmethod), full.names = TRUE)
datfile <- datfiles[length(datfiles)]
print(paste("Loading this data file: ", datfile))
load(datfile)
load(file.path(datdir, pasteu(exptstr, "se_filtered.Rda")))



if (opt$clusmethod == "rsec"){
  se_filtered <- se_filtered[, colnames(cl)]
  dat <- log2(assay(cl) + 1)
  colData <- DataFrame(cl@clusterMatrix, expt = colData(se_filtered)$expt,
                                      batch = colData(se_filtered)$batch,
                      samples = colnames(cl))
} else if (opt$clusmethod == "snn"){
  se_filtered <- se_filtered[, colnames(seu@data)]
  dat <- seu@data
  metadata <- data.frame(seu@meta.data, expt = colData(se_filtered)$expt, 
	    			      batch = colData(se_filtered)$batch,  
	    	      samples = colnames(seu@data))
metadata <- metadata[with(metadata, order(res.2, expt, batch)),]
#metadata <- metadata[,c(grep('^res|expt|batch', colnames(metadata)))]
metadata <- metadata[,c('res.2', 'res.0.5', 'expt', 'batch')]

}

seed <- 2782472

# heatmaps
bigPalette <- bigPalette[-3]
colRKC<- c("chartreuse3", "firebrick3", "chocolate4","slategray2","darkviolet", 
                "darkorange2","pink", "gold", "deepskyblue", "pink3", 
                "deeppink", "grey36", "royalblue3", "mediumorchid1","grey0", 
                "cyan2", "darkseagreen4")
clusIndex <- c('1' ,'2', '3','4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16')
#clusterLegend <- setNames(rep(list(bigPalette), ncol(metadata)), colnames(metadata))
clusterLegend <- setNames(rep(list(colRKC), ncol(colData)),colnames(colData))
#orderedLegend  <- factor(mydata(metadata)[,"mergeClusters"],
#                      levels=clusterLegend(metadata)[["mergeClusters"]][, "name"])
#names(orderedLegend) <- colnames(metadata)
clusterLegend[["expt"]] <- cole

markers <- intersect(unlist(read.table(file.path("../ref", opt$markerfile))), 
                        rownames(dat))
breakv <- c(min(dat), 
            seq(0, quantile(dat[dat > 0], .99, na.rm = TRUE), length = 50), 
            max(dat))
breakv <- unique(breakv)
pdf(file = file.path(vizdir, pasteu0(exptstr, "markerhm", method, opt$norm, opt$clusmethod, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")), 
    width = 8.5, height = 11)
if (opt$clusmethod == "snn") {
   plotHeatmap(dat[markers, rownames(metadata)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = metadata, clusterSamplesData = "dendrogramValue",
				 clusterLegend = clusterLegend, 
				 annLegend = TRUE, overRideClusterLimit = TRUE)
} else if (opt$clusmethod == "rsec"){
  plotHeatmap(cl, clusterSamples = TRUE, clusterFeaturesData = markers, clusterFeatures = FALSE, breaks = breakv, 
				 clusterLegend = list(expt = cole), colData = colData,
				 annLegend = TRUE, overRideClusterLimit = TRUE)

}
dev.off()
###########################################
# t-SNE colored by cluster and time point

if (opt$clusmethod == "snn"){
#if (is.null(seu@dr$tsne)){
print("recomputing t-SNE")
#ngenesvec = c(500, 1000, 3000, 5000, nrow(seu@data))
ngenesvec = c(500)
#perpvec = c(50, 60, 70, 80, 90)
perpvec= c(80)
vars <- apply(dat, 1, var)
vars <- sort(vars, decreasing = TRUE)

for (ngenes in ngenesvec){
for (perp in perpvec){
genes.use <- names(vars)[1:ngenes]

if (method == "scone"){
  seu <- RunTSNE(seu, reduction.use = "pca", dims.use = 1:50,
  genes.use = genes.use, seed.use = seed, tsne.method = "Rtsne", perplexity = perp, max_iter = 10000,
  dim.embed = 2, reduction.name = "tsne")
} else if (method == "zinb"){
  seu <- RunTSNE(seu, reduction.use = "zinbwave", dims.use = 1:20,
  genes.use = genes.use, seed.use = seed, tsne.method = "Rtsne", perplexity = perp, max_iter = 10000,
  dim.embed = 2, reduction.name = "tsne")
}
#save(seu, file = datfile)

pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", ngenes, "p", perp, "clus", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
plot(seu@dr$tsne@cell.embeddings, pch = 19, cex = 0.4, col = alpha(colRKC[factor(seu@meta.data[,"res.2"])], 0.3), 
				  xlab = "TSNE 1", ylab = "TSNE 2", main = paste(ngenes, "genes, perplexity =", perp))
legend("topleft", legend = levels(factor(seu@meta.data[,"res.2"])), fill = colRKC, cex = 0.5)
dev.off()
}
}
#}
}

# if (opt$clusmethod == "rsec"){ 
# if (is.null(cl@reducedDims$tsne)){
# ngenes = nrow(cl)
# vars <- apply(dat, 1, var)
# vars <- sort(vars, decreasing = TRUE)
# genes.use <- names(vars)[1:ngenes]

# if (method == "scone"){
#   seu <- RunTSNE(seu, reduction.use = "pca", dims.use = 1:50,
#   genes.use = genes.use, seed.use = seed, tsne.method = "Rtsne", perplexity = 10, max_iter = 10000,
#   dim.embed = 2, reduction.name = "tsne")
# } else if (method == "zinb"){
#   seu <- RunTSNE(seu, reduction.use = "zinbwave", dims.use = 1:20,
#   genes.use = genes.use, seed.use = seed, tsne.method = "Rtsne", perplexity = 10, max_iter = 10000,
#   dim.embed = 2, reduction.name = "tsne")
# }
# save(seu, file = datfile)
# }}

#  pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", "clus", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
#  plot(seu@dr$tsne@cell.embeddings, pch = 19, cex = 0.4, col = alpha(colRKC[factor(seu@meta.data[,"res.2"])], 0.3), 
#  				  xlab = "TSNE 1", ylab = "TSNE 2")
#  legend("topleft", legend = levels(factor(seu@meta.data[,"res.2"])), fill = colRKC, cex = 0.5)
#  dev.off()

# pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", "batch", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
# plot(seu@dr$tsne@cell.embeddings, pch = 19, cex = 0.4, col = alpha(colb[batch], 0.3), xlab = "TSNE 1", ylab =" TSNE 2")
# legend("bottomleft", legend = levels(batch), fill = colb, cex = 0.6)
# dev.off()

# pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", "expt", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
# plot(seu@dr$tsne@cell.embeddings, pch = 19, cex = 0.4, col = alpha(cole[expt], 0.3), xlab = "TSNE 1", ylab =" TSNE 2")
# legend("bottomleft", legend = levels(expt), fill = cole, cex = 0.6)
# dev.off()

#  pdf(file = file.path(vizdir, pasteu0(exptstr, "tsne", "geneexp", method, opt$norm, format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")))
#  t1 <- theme(plot.background=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(),axis.ticks=element_blank(), legend.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),legend.key= element_rect(fill="white"), panel.border = element_rect(fill=NA,colour = "black"),axis.line=element_blank(),aspect.ratio=1)
#  dat <- data.frame(seu@dr$tsne@cell.embeddings, t(seu@data[markers,]))
#  par(mar=c(2,2,1,1), mfrow=c(1,1))
#  for (gene in markers){
#    p <- ggplot(dat, aes_string("tSNE_1", "tSNE_2", colour = gene)) + geom_point(cex=0.5)
#    print(p + t1 +
#            scale_colour_gradient2(low = "#053061", mid = "grey95", high = "#67001F") +
#            ggtitle(gene))
#  }
#  dev.off()