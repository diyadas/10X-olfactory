# Filtering of 10X data
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: June 1, 2018

# Load command-line arguments
rm(list=ls()); options(getClass.msg=FALSE)
library(optparse)
option_list <- list(
  make_option("--expt", type="character", help="Experiment ID"),
  make_option("--ncores", default="1", type="double"),
  make_option("--aggr", type="character", help="name of aggregate ID/directory"),
  make_option("--exptinfo", type="character", help="text file with list of cellranger experiments, from aggregate...csv, expanded with columns for expt and batch: \n column 1: library_id, column 2: expt, column 3: batch"),
  make_option("--annotation", default="GRCm38p4Mm10", type=character, help="name of genomic reference used by cellranger")
)
opt <- parse_args(OptionParser(option_list=option_list))
exptstr <- opt$expt
outdir <- file.path("../output", expt, "data")
vizdir <- file.path("../output", expt, "viz")
crdir <- file.path("../output", expt, "cr")
exptinfo <- read.csv(file.path("../output", expt, opt$exptinfo), stringsAsFactors = FALSE)

# Set up packages and parallel environment
library(BiocParallel)
register(MulticoreParam(workers = opt$ncores))
library(cellrangerRkit)
library(SummarizedExperiment)

library(scater)
library(scone)
library(clusterExperiment)
library(ggplot2)
library(magrittr)
library(cowplot)
library(Rtsne)
library(zinbwave)

library(scales)

# Custom code
# Create color palettes
pal <- clusterExperiment::bigPalette
colb <- c("dodgerblue3", "darkviolet", "goldenrod", "chartreuse3", "cadetblue2", "magenta", "chocolate1", "forestgreen", "plum", "azure3", "chocolate4", "darkslategray", "cornsilk", "aquamarine3", "burlywood3", "darkblue", "gray45")
cole <- c("cornflowerblue", "darkgoldenrod", "darkorchid", "darkorange2", "deeppink3", "cadetblue1", "azure4", "darkslateblue", "darkolivegreen1", "antiquewhite2")
fastpca <- function(expr, scale=FALSE) {
  k <- 50
  svd_raw <- rARPACK::svds(scale(t(expr), center=TRUE, scale=scale), k=k, nu=k, nv=0)
  pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:k])
  return(pc_raw)
}


# ---- Loading data ----
# Load counts, batch, expt and create SummarizedExperiment
expt <- as.factor(as.vector(apply(exptinfo, 1, function(crexpt){
  crmat <- load_cellranger_matrix(file.path(crdir,crexpt[1])) # crexpt[1] is the library_id field, i.e. RCOB2A
  return(rep(crexpt[2], ncol(exprs(crmat))))
})))

batch <- as.factor(as.vector(apply(exptinfo, 1, function(crexpt){
  crmat <- load_cellranger_matrix(file.path(crdir,crexpt[1])) # crexpt[1] is the library_id field, i.e. RCOB2A
  return(rep(crexpt[3], ncol(exprs(crmat))))
})))
  
counts <- as.matrix(exprs(load_cellranger_matrix(file.path(crdir,opt$aggr))))

se <- SummarizedExperiment(list(counts=counts), colData=data.frame(batch=batch, expt=expt))
genes <- read.table(file=file.path(crdir, opt$aggr, "/outs/filtered_gene_bc_matrices_mex",opt$annotation,"genes.tsv"))
colnames(genes) <- c("geneID", "Symbol")
rowData(se) <- genes
head(rowData(se))
head(colData(se))


# Exploratory Data Analysis

hist(colSums(assay(se)), breaks=30, 
     xlab='Number of UMI', main="Number of UMI per sample")

hist(colSums(assay(se)>0)/nrow(se), breaks=30, 
     xlab="Proportion of detected genes", main="Proportion of detected genes")

hist(colSums(assay(se)>0), breaks=30, 
     xlab="Number of detected genes", main="Number of detected genes")

hist(colMeans(assay(se)==0), breaks=30, 
     xlab="Proportion of zeros", main="Proportion of zeros")

boxplot(colSums(assay(se))~colData(se)$batch, main="Number of UMI per sample", col=colb, las=2)
boxplot(colSums(assay(se)>0)~colData(se)$batch, main="Number of detected genes", col=colb, las=2)
boxplot(colMeans(assay(se)==0)~colData(se)$batch, main="Proportion of zeros", col=colb, las=2)

simple <- assay(se)[rowSums(assay(se))>10,]
simple <- SUM_FN(simple)

runQCmetrics <- FALSE
if (runQCmetrics) {
  pca <- prcomp(t(simple), scale. = TRUE)
  # pca <- fastpca(log2(simple + 0.1))
  ###---Use scater package to calculate some quality control metrics. Assess QC metrics.
  
  ###only works with unique identifiers - Ensembl ids, not gene names b/c some gene name duplication (65/28000)
  sce <- newSCESet(countData = assay(se))
  sce <- calculateQCMetrics(sce)
  
  rownames(se) <- rowData(se)$Symbol
  head(rownames(se))
  
  ###only works with gene names
  ribo_idx <- grep("^Rpl|^Rps", rowData(se)[,2])
  mito_idx <- grep("^Mt", rowData(se)[,2])
  ribo_pct <- colSums(assay(se)[ribo_idx,])/colSums(assay(se)) * 100
  mito_pct <- colSums(assay(se)[mito_idx,])/colSums(assay(se)) * 100
  plot(mito_pct, col=colb[batch], xlab="cell index", main="% mito (Mt*) genes"); legend("topleft", legend=levels(batch), fill=colb, cex=0.8)
  boxplot(mito_pct ~ colData(se)$batch, main="percent mito genes", col=colb, las=2, cex.axis=0.7)
  plot(ribo_pct, col=colb[batch], xlab="cell index", main="% ribo (Rpl*) genes"); legend("topleft", legend=levels(batch), fill=colb, cex=0.8)
  boxplot(ribo_pct ~ colData(se)$batch, main="percent ribo genes", col=colb, las=2, cex.axis=0.7)
  
  qc <- as.matrix(data.frame(colData(sce)[,c(2, 4:8)], mito_pct = mito_pct, ribo_pct = ribo_pct))
  
  qcpca <- prcomp(qc, scale. = TRUE)
  screeplot(qcpca, type = "lines", main = "QC-PCA screeplot")
  plot(qcpca$x, col=alpha(colb[batch], 0.3), pch=19, main="QC PCA")
  legend("topleft", legend=levels(batch), fill=colb, cex=0.5)
  
  rm(mito_idx, ribo_idx, fig_pca, fig_qpca)
  save(se, sce, qc, pca, qcpca, file=paste0("dataObjects/", exptStr, "_unfiltered.Rda"))
} else {
  load(paste0("dataObjects/", exptStr, "_unfiltered.Rda"))
}

boxplot(qc[,7] ~ colData(se)$batch, main="percent mito genes", col=colb, las=2, cex.axis=0.7)
boxplot(qc[,8] ~ colData(se)$batch, main="percent ribo genes", col=colb, las=2, cex.axis=0.7)

screeplot(qcpca, type = "lines", main = "QC-PCA screeplot")
plot(qcpca$x, col=alpha(colb[batch], 0.3), pch=19, main="QC PCA")
legend("topleft", legend=levels(batch), fill=colb, cex=0.5)

screeplot(pca, type = "lines", npcs = 50, main = "Expression PCA screeplot")
plot(pca$x, col=alpha(colb[batch], 0.3), pch=19, cex=0.4, main="simple filtering, scaling --> PCA")
legend("topright", legend=levels(batch), fill=colb, cex=0.5)

print(paste0("Percent total variance captured in first 10 expression PCs is ", round(cumsum(pca$sdev^2/sum(pca$sdev^2))[10] * 100, digits = 2)))

fig_data <- data.frame(pca$x[,1:2], qc,
                       QPC1=qcpca$x[,1], QPC2=qcpca$x[,2])

fig_pca <- ggplot(fig_data, aes(x = PC1, y = PC2, color = log10_total_counts)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_pca

fig_qpca <- ggplot(fig_data, aes(x = QPC1, y = QPC2, color = log10_total_counts)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_qpca

fig_pca <- ggplot(fig_data, aes(x = PC1, y = PC2, color = ribo_pct)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_pca

fig_qpca <- ggplot(fig_data, aes(x = QPC1, y = QPC2, color = ribo_pct)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_qpca

fig_pca <- ggplot(fig_data, aes(x = PC1, y = PC2, color = mito_pct)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_pca

fig_qpca <- ggplot(fig_data, aes(x = QPC1, y = QPC2, color = mito_pct)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_qpca

fig_pca <- ggplot(fig_data, aes(x = PC1, y = PC2, color = log10_total_features)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_pca

fig_qpca <- ggplot(fig_data, aes(x = QPC1, y = QPC2, color = log10_total_features)) +
  geom_point() + scale_color_continuous(low = "blue", high = "yellow")
fig_qpca



# Filtering

colData(se) <- cbind(colData(se), qc)

# select expressed genes only
dim(se)
se <- se[rowSums(assay(se))>0,]
dim(se)

# select common genes
num_reads <- quantile(assay(se)[assay(se) > 0])[4]
num_cells = 0.25*ncol(se)
is_common = rowSums(assay(se) >= num_reads) >= num_cells
table(is_common)

#hk <- read.table("ref/HK_genes.txt")
hk <- read.table("ref/hkl615.txt")
hk <- as.character(hk[,1])
# hk <- sample(hk, 1000)
# write.table(hk, file="ref/HK_genes_1000sampled.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
hk <- intersect(hk, rowData(se)$Symbol)
hk_idx <- which(rowData(se)$Symbol %in% hk)

mfilt <- metric_sample_filter(assay(se),
                              nreads = colData(sce)$total_counts,
                              gene_filter = is_common,
                              pos_controls = hk_idx,
                              hard_nreads = 2000,
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)
pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_mfiltOutput_zcut3_hnr2K.pdf"), title = "metric_sample_filtering")
mfilt <- metric_sample_filter(assay(se),
                              nreads = colData(sce)$total_counts,
                              gene_filter = is_common,
                              pos_controls = hk_idx,
                              hard_nreads = 2000,
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)
dev.off()

plot(pca$x, pch=19, col=pal[as.numeric(mfilt$filtered_breadth)+1],
     main = "PCA Filtered on transcriptome 'breadth'")
plot(pca$x, pch=19, col=pal[as.numeric(mfilt$filtered_fnr)+1],
     main = "PCA Filtered on FNR AUC")
plot(pca$x, pch=19, col=pal[as.numeric(mfilt$filtered_nreads)+1],
     main = "PCA Filtered on nreads")

plot(qcpca$x, pch=19, col=pal[as.numeric(mfilt$filtered_breadth)+1],
     main = "QPCA Filtered on transcriptome 'breadth'")
plot(qcpca$x, pch=19, col=pal[as.numeric(mfilt$filtered_fnr)+1],
     main = "QPCA Filtered on FNR AUC")
plot(qcpca$x, pch=19, col=pal[as.numeric(mfilt$filtered_nreads)+1],
     main = "QPCA Filtered on nreads")

table(mfilt$filtered_nreads, mfilt$filtered_fnr)
filter_cell <- !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)

plot(qcpca$x, pch=19, col=pal[as.numeric(filter_cell)+1],
     main = "PCA Filtered")


# Final Gene Filtering: Highly expressed in at least 5 cells
num_reads <- quantile(assay(se)[assay(se) > 0])[4]
num_cells = 3
is_quality = rowSums(assay(se) >= num_reads ) >= num_cells
table(is_quality)
filtered <- se[is_quality, filter_cell]
dim(filtered)
rownames(filtered) <- rowData(filtered)$Symbol
qc <- qc[colnames(filtered),]
batch <- colData(filtered)$batch
expt <- colData(filtered)$expt

# Check for batch effects prior to normalization

# check for batch effects
plot(pca$x, pch=19, col=colb[colData(se)$batch],
     main = "PCA Color-coded by batch")
legend("topleft", legend=levels(batch), fill=colb, cex=0.6)
plot(qcpca$x, pch=19, col=colb[colData(se)$batch],
     main = "QPCA Color-coded by batch")
legend("bottomleft", legend=levels(batch), fill=colb, cex=0.6)

boxplot(pca$x[,1] ~ colData(se)$batch, col=colb, cex.axis = 0.75, las = 2, main = "First Principal Component")
boxplot(pca$x[,2] ~ colData(se)$batch, col=colb, cex.axis = 0.75, las = 2, main = "Second Principal Component")

simple <- assay(filtered)

pca <- fastpca(log2(simple + 1))
tsne_data <- Rtsne(pca[,1:10], pca=FALSE, max_iter = 5000)
plot(tsne_data$Y, pch=19, cex=0.4, col=alpha(colb[colData(filtered)$batch],0.5))
legend("topleft", legend=levels(batch), fill=colb, cex=0.5)
pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_tSNE_50PCs_filteredPreNorm.pdf"), title = "tSNE_50PCs_filteredPreNorm")
plot(tsne_data$Y, pch=19, cex=0.4, col=alpha(colb[colData(filtered)$batch],0.5))
legend("topleft", legend=levels(batch), fill=colb, cex=0.5)
dev.off()


# Preparation for normalization
oeDEgenes <- read.table(file="ref/oeRegPosCon.txt")
poscon <- intersect(oeDEgenes$V1, rowData(filtered)$Symbol)
rowData(filtered)$poscon <- (rowData(filtered)$Symbol %in% poscon)

## select negative controls (housekeeping)
hk <- intersect(hk, rowData(filtered)$Symbol)

negconeval <- sample(hk, length(poscon)/2)
negconruv <- setdiff(hk, negconeval)

rowData(filtered)$negcon_eval <- (rowData(filtered)$Symbol %in% negconeval)
rowData(filtered)$negcon_ruv <- (rowData(filtered)$Symbol %in% negconruv)



# Gene expression plots post filtering, but pre-normalization

filtCounts <- assay(filtered)
logfiltCounts <- log2(filtCounts+1)

pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_preNorm_poscon_heatmap.pdf"), title = "poscon heatmap")
plotHeatmap(logfiltCounts[rowData(filtered)$poscon,], sampleData=data.frame(expt = colData(filtered)[,2], batch = colData(filtered)[,1]), clusterLegend = list(expt = cole, batch=colb), main="Positive controls", breaks=.99)
# aheatmap(logfiltCounts[rowData(filtered)$poscon,], annCol=data.frame(batch = colData(filtered)[,1]),  main="Positive controls", breaks=.99)
dev.off()
plotHeatmap(logfiltCounts[rowData(filtered)$poscon,], sampleData=data.frame(batch = colData(filtered)[,1], expt = colData(filtered)[,2]), clusterLegend = list(batch=colb, expt = cole), main="Positive controls", breaks=.99)

pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_preNorm_negconeval_heatmap.pdf"), title = "negconeval heatmap")
plotHeatmap(logfiltCounts[rowData(filtered)$negcon_eval,], sampleData=data.frame(batch = colData(filtered)[,1], expt = colData(filtered)[,2]), clusterLegend = list(batch=colb, expt = cole), main="Negative controls", breaks=.99)
dev.off()
plotHeatmap(logfiltCounts[rowData(filtered)$negcon_eval,], sampleData=data.frame(batch = colData(filtered)[,1], expt = colData(filtered)[,2]), clusterLegend = list(batch=colb, expt = cole),  main="Negative controls", breaks=.99)

pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_preNorm_negconruv_heatmap.pdf"), title = "negconruv heatmap")
plotHeatmap(logfiltCounts[rowData(filtered)$negcon_ruv,], sampleData=data.frame(batch = colData(filtered)[,1], expt = colData(filtered)[,2]), clusterLegend = list(batch=colb, expt = cole), main="Negative controls", breaks=.99)
dev.off()

##Correlation of QC metrics with expression PCs for filtered, pre-norm 
cors <- lapply(1:5, function(i) abs(cor(pca[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=factor(rep(colnames(qc), 5), levels=colnames(qc)),
                   Dimension=as.factor(rep(paste0("PC", 1:5), each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=bigPalette) + ylim(0, 1) + 
  ggtitle("Correlation between QC and expression PCA") -> fig_barplot
fig_barplot
pdf(file=paste0("output/viz/filteredPreNorm/", exptStr, "_filteredPreNorm_corQCandExprPCA.pdf"))
fig_barplot
dev.off()
###################

save(qc, pca, filtered, batch, expt, file=paste0("dataObjects/", exptStr, "_filtered.rda"))
save(filtCounts, logfiltCounts, file=paste0("dataObjects/", exptStr, "_filteredCounts.rda"))



# SCONE Normalization
## preparing to run SCONE to test/rank normalizations
#####run SCONE on the cluster

run_scone <- TRUE
if(run_scone) {
scone_obj <- SconeExperiment(filtered,
               which_qc = which(colnames(colData(filtered)) %in% colnames(qc)),
               which_negconruv = 5L,
               which_negconeval = 4L,
               which_poscon = 3L,
               which_batch = 1L,
               which_bio = 2L)

save(scone_obj, file=paste0("dataObjects/", exptStr, "_scone_obj1.Rda"))
} else { 
  load(paste0("dataObjects/", exptStr, "_scone_normMats.Rda"))
}

#SCONE Normalization
###SCONE output evaluation, choosing normalizations

eval_scone <- TRUE
if(eval_scone) {
load(paste0("~/Research/Projects/SCRNASeq/10XG/regen/dataObjects/", exptStr, "_scone_res_sub10.Rda")
sconeRanks <- sort(rowMeans(sapply(sub_ranks, function(x) x[[2]][names(sub_ranks[[1]][[2]])])), decreasing = TRUE)
################
################

sconeObjParams <- scone(scone_obj, scaling = list(none=identity, sum = SUM_FN, tmm = TMM_FN,
                                   uq = UQ_FN,fq = FQT_FN, deseq = DESEQ_FN), k_ruv = 3,
      k_qc = 5, adjust_bio = "no", adjust_batch = "yes", eval_kclust = 5:15, zero = "postadjust", run=FALSE)

# scoreRanks <- get_score_ranks(scone_obj2)
scone_obj2 <- scone_obj
scone_obj2@scone_params <- sconeObjParams@scone_params[names(sconeRanks)[1:6],] 

save(scone_obj2, file=paste0(dataDir, exptStr, "_scone_obj2.Rda"))
} 

#####-----running scone on selected normalization: none,fq,qc_k=2,no_bio,batch and retrieving it in memory
#####
#load("regen2_scone_normMats.Rda")



###Prep/Run ZINB-WaVE using most variable genes based on filtered se

library(Seurat)
library(Rtsne)
library(zinbwave)
library(clusterExperiment)
#library(cellrangerRkit)

#NMF::nmf.options(grid.patch=TRUE)##put anywhere in the script to remove first (blank) page of heatmap pdfs
importFiltered <- FALSE
if (importFiltered){
#pal <- clusterExperiment::bigPalette
load(paste0("dataObjects/", exptStr, "_filtered.rda"))
load(paste0("dataObjects/", exptStr, "_filteredCounts.rda"))
# rownames(filtered) <- rowData(filtered)$Symbol
# batch <- colData(filtered)$batch
}


# ZINB-WaVE set-up

# filtCounts <- assay(filtered)
# logfiltCounts <- log2(filtCounts+1)

# vars <- rowVars(logSconeNorm)
# names(vars) <- rownames(logSconeNorm)
vars <- rowVars(logfiltCounts)
names(vars) <- rownames(logfiltCounts)
vars <- sort(vars, decreasing = TRUE)
vargenes <- names(vars)[1:1000]

save(vargenes, filtered, file=paste0("dataObjects/", exptStr, "_forZinbwaveRun.Rda"))

###running ZINB-WaVE and visualizing output -- running ZINB-WaVE on cluster 
run_zinbwave <- FALSE
if(run_zinbwave) {
library(BiocParallel)
library(doParallel)
registerDoParallel(3)
register(DoparParam())

zinb <- zinbFit(filtered[vargenes,], X = "~ log10_total_counts + ribo_pct", K=10, epsilon=1000)
save(zinb, file=paste0("dataObjects/", exptStr, "_zinbwave2Kvar_TC_ribo.rda"))
W <- getW(zinb)
save(W, file = paste0("dataObjects/", exptStr, "_zinbW.Rda"))
d <- dist(W)
tsne_zinb <- Rtsne(d, is_distance = TRUE, pca = FALSE, max_iter=5000)
} else {
  # load("dataObjects/regen2_zinbwave2Kvar_TC_ribo.Rda")
  # load("dataObjects/regen2_zinbW.Rda")
  # load("dataObjects/regen2_tsne_zinb.Rda")
  load(paste0("dataObjects/", exptStr, "_zinbwave1Kvar_TC_top500_ribo.Rda"))
  load(paste0("dataObjects/", exptStr, "_zinbW.Rda"))
  load(paste0("dataObjects/", exptStr, "_tsne_zinb.Rda"))
}

plot(tsne_zinb$Y, pch=19, cex=0.4, col=colb[batch], main = "t-SNE (ZINB-WaVE), 1KvarGenes")
legend("topleft", legend=levels(batch), fill=colb, cex=0.35)
pdf(file = paste0("output/viz/", exptStr, "_zinb_tSNE_batch.pdf"))
plot(tsne_zinb$Y, pch=19, cex=0.4, col=colb[batch], main = "t-SNE (ZINB-WaVE), 1KvarGenes")
legend("topleft", legend=levels(batch), fill=colb, cex=0.35)
dev.off()
plot(tsne_zinb$Y, pch=19, cex=0.4, col=alpha(cole[expt], 0.25), main = "t-SNE (ZINB-WaVE), 1KvarGenes")
legend("topleft", legend=levels(expt), fill=cole, cex=0.5)
pdf(file = paste0("output/viz/", exptStr, "_zinb_tSNE_expt.pdf"))
plot(tsne_zinb$Y, pch=19, cex=0.4, col=alpha(cole[expt], 0.25), main = "t-SNE (ZINB-WaVE), 1KvarGenes")
legend("topleft", legend=levels(expt), fill=cole, cex=0.5)
dev.off()
# plot(W, pch=19, col=pal[batch])
# pairs(W[,1:5], pch=19, col=pal[batch])


###ZINB-WaVE QC correlation
#####Correlation of QC metrics with ZINB-WaVE W projection components

zbrd <- getW(zinb)
cors <- lapply(1:5, function(i) abs(cor(zbrd[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=factor(rep(colnames(qc), 5), levels=colnames(qc)),
                   Dimension=as.factor(rep(paste0("zinbW", 1:5), each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=bigPalette) + ylim(0, 1) + 
  ggtitle("Correlation between QC and ZINBWaVE W (2KVar)") -> fig_barplot
fig_barplot
pdf(file="output/viz/zinb2KTCtop50Ribo_corQCandZinbW_.pdf")
fig_barplot
dev.off()

##Evaluating SCONE subSampling output

load(paste0(dataDir, exptStr, "_scone_res_sub10.rda"))
sconeRanks <- sort(rowMeans(sapply(sub_ranks, function(x) x[[2]][names(sub_ranks[[1]][[2]])])),
     decreasing = TRUE)



##SCONE Normalization
###Retrieving chosen normalizations

load("dataObjects/regen2_scone_normMats.Rda")
sconeNorm <- get_normalized(scone_obj4, "none,fq,qc_k=2,no_bio,no_batch")
rownames(sconeNorm) <- rowData(filtered)[,2]
logSconeNorm <- log2(sconeNorm+1)
save(sconeNorm, logSconeNorm, file="dataObjects/regen2_normCounts.Rda")

#######
#load("dataObjects/regen2_scone_normMats.Rda")
sconeNorm2 <- get_normalized(scone_obj4, "none,fq,qc_k=3,no_bio,no_batch")
rownames(sconeNorm2) <- rowData(filtered)[,2]
logSconeNorm2 <- log2(sconeNorm2+1)
save(sconeNorm2, logSconeNorm2, file="dataObjects/regen2_normCounts2.Rda")

# #######
# load("dataObjects/regen2_scone_normMats.Rda")
sconeNorm3 <- get_normalized(scone_obj4, "none,fq,qc_k=2,no_bio,batch")
rownames(sconeNorm3) <- rowData(filtered)[,2]
logSconeNorm3 <- log2(sconeNorm3+1)
save(sconeNorm3, logSconeNorm3, file="dataObjects/regen2_normCounts3.Rda")

# #######
# #load("dataObjects/regen2_scone_normMats.Rda")
# sconeNorm4 <- get_normalized(scone_obj4, "none,fq,qc_k=3,no_bio,batch")
# rownames(sconeNorm4) <- rowData(filtered)[,2]
# logSconeNorm4 <- log2(sconeNorm4+1)
#save(sconeNorm4, logSconeNorm4, file="dataObjects/regen2_normCounts4.Rda")

# biplot_interactive(scone_obj)


##Correlation of QC metrics with expression PCs for scone norm 

# pca_sconeNorm <- prcomp(t(logSconeNorm))
pca_sconeNorm <- fastpca(logSconeNorm)

cors <- lapply(1:5, function(i) abs(cor(pca_sconeNorm[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=factor(rep(colnames(qc), 5), levels=colnames(qc)),
                   Dimension=as.factor(rep(paste0("PC", 1:5), each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=bigPalette) + ylim(0, 1) + 
  ggtitle("Correlation between QC and expression PCA") -> fig_barplot
fig_barplot
pdf(file="output/viz/sconeNorm_allGenes_corQCandExprPCA.pdf")
fig_barplot
dev.off()
###################

##Correlation of QC metrics with expression PCs for additional scone norms 
# pca_sconeNorm2 <- prcomp(t(logSconeNorm2))
pca_sconeNorm2 <- fastpca(logSconeNorm2)

cors <- lapply(1:5, function(i) abs(cor(pca_sconeNorm2[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=factor(rep(colnames(qc), 5), levels=colnames(qc)),
                   Dimension=as.factor(rep(paste0("PC", 1:5), each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=bigPalette) + ylim(0, 1) + 
  ggtitle("Correlation between QC and expression PCA") -> fig_barplot
fig_barplot
pdf(file="output/viz/sconeNorm2_allGenes_corQCandExprPCA.pdf")
fig_barplot
dev.off()
###################

# pca_sconeNorm3 <- prcomp(t(logSconeNorm3))
pca_sconeNorm3 <- fastpca(logSconeNorm3)

cors <- lapply(1:5, function(i) abs(cor(pca_sconeNorm3[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=factor(rep(colnames(qc), 5), levels=colnames(qc)),
                   Dimension=as.factor(rep(paste0("PC", 1:5), each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=bigPalette) + ylim(0, 1) + 
  ggtitle("Correlation between QC and expression PCA") -> fig_barplot
fig_barplot
pdf(file="output/viz/sconeNorm3_allGenes_corQCandExprPCA.pdf")
fig_barplot
dev.off()
