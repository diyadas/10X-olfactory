# Filtering of 10X data
# Authors: Russell Fletcher, Diya Das, and Rebecca Chance
# Last revised: Tue Jun 19 11:54:38 2018

# Load command-line arguments
rm(list = ls())
options(getClass.msg = FALSE)
library(optparse)

option_list <- list(
  make_option("--expt", type = "character", help = "Experiment ID"),
  make_option("--ncores", default = "1", type = "double"),
  make_option("--aggr", type = "character",
              help = "name of aggregate ID/directory"),
  make_option("--exptinfo", type = "character",
              help = "text file with list of cellranger experiments, from 
              aggregate...csv, expanded with columns for expt and batch: \n 
              column 1: library_id, column 3: expt, column 4: batch"),
  make_option("--annotation", default = "GRCm38p4Mm10", type = "character",
              help = "name of genomic reference used by cellranger"),
  make_option("--hkfile", default = "../ref/hkl615.txt", type = "character", 
              help = "path to file containing housekeeping genes"),
  make_option("--posctrlfile", default = "../ref/oeRegPosCon.txt", 
              type = "character",
              help = "path to file containing positive control genes"),
  make_option("--runQC", default = FALSE, type = "logical",
  	      help = "whether to calculate QC metrics, will error if
	      QC metrics have not previously been calculated"),
  make_option("--fast", default = FALSE, type = "logical",
  	      help = "whether to use fast (approximate) PCA")
  )
opt <- parse_args(OptionParser(option_list = option_list))
exptstr <- opt$expt
outdir <- file.path("../output", exptstr, "data")
vizdir <- file.path("../output", exptstr, "viz", "prenormalization")
crdir <- file.path("../output", exptstr, "cr")
exptinfo <- read.csv(file.path("../output", exptstr, opt$exptinfo),
                     stringsAsFactors = FALSE)
system(paste("mkdir -p", c(outdir, vizdir)))

# Set up packages and parallel environment
library(BiocParallel)
register(MulticoreParam(workers = opt$ncores))
library(cellrangerRkit)
library(SummarizedExperiment)
library(scone)
library(scater)
library(Rtsne)
library(ggplot2)
library(scales)
library(clusterExperiment)

# Source helper functions
source("tenx_helper.R")

# ---- Loading data ----
# Load counts, batch, expt and create SummarizedExperiment
expt_batch <- apply(exptinfo, 1, function(crexpt) {
  crmat <- load_cellranger_matrix(file.path(crdir, crexpt[1]))
  # crexpt[1] is the library_id field, i.e. RCOB2A
  return(data.frame(expt = rep(crexpt[3], ncol(exprs(crmat))), 
                    batch = rep(crexpt[4], ncol(exprs(crmat)))))
})
expt_batch <- do.call(rbind, expt_batch)
expt <- expt_batch$expt
batch <- expt_batch$expt

counts <- as.matrix(exprs(load_cellranger_matrix(file.path(crdir, opt$aggr))))

se <- SummarizedExperiment(list(counts = counts), 
                           colData = data.frame(batch = batch, expt = expt))
genes <- read.table(file = file.path(crdir, opt$aggr, 
                                     "outs/filtered_gene_bc_matrices_mex",
                                     opt$annotation, "genes.tsv"))
colnames(genes) <- c("geneID", "Symbol")
rowData(se) <- genes
rownames(se) <- rowData(se)$Symbol

# Exploratory Data Analysis
pdf(file = file.path(vizdir, pasteu(exptstr, "EDA_prefilt.pdf")),
    height = 11, width = 8.5)
detectedgeneplot(se, "Pre-filtering\n") # function in helper script
dev.off()

# ---- Calculate QC metrics and plot ----
runQC <- opt$runQC
fast <- opt$fast
if (runQC) {
  se_simple <- assay(se)[rowSums(assay(se)) > 10, ]
  se_simple <- SUM_FN(se_simple)
  if (fast){
     pca <- fastpca(log2(se_simple + 0.1))
  } else {
     pca <- prcomp(t(log2(se_simple + 0.1)))
  }
  # Use scater package to calculate some quality control metrics.
  # only works with unique identifiers - Ensembl IDs,
  # not gene names b/c some gene name duplication (~65/28000)
  sce <- as(se, "SingleCellExperiment")
  rownames(sce) <- rowData(sce)$geneID
  sce <- calculateQCMetrics(sce)
  
  # Calculate mitochondrial and ribosomal percentage in each sample
  ribo_idx <- grep("^Rpl|^Rps", rowData(se)$Symbol)
  mito_idx <- grep("^Mt", rowData(se)$Symbol)
  ribo_pct <- colSums(assay(se)[ribo_idx, ]) / colSums(assay(se)) * 100
  mito_pct <- colSums(assay(se)[mito_idx, ]) / colSums(assay(se)) * 100
  
  qc <- as.matrix(data.frame(colData(sce)[, c("total_features_by_counts", 
     					      "total_counts",
					      "log10_total_counts",
					      "pct_counts_in_top_50_features",
					      "pct_counts_in_top_100_features",
					      "pct_counts_in_top_200_features",
					      "log10_total_features")], 
			     mito_pct = mito_pct, 
                  	     ribo_pct = ribo_pct))
  qcpca <- prcomp(qc, scale. = TRUE)
 
 pdf(file = file.path(vizdir, pasteu(exptstr, "mitoribo_prefilt.pdf")),
    height = 11, width = 8.5)
plot(qc[,"mito_pct"], col = colb[batch], xlab = "cell index",
     main = "% mito (Mt*) genes")
legend("topleft", legend = levels(batch), fill = colb, cex = 0.8)
boxplot(qc[,"mito_pct"] ~ colData(se)$batch, main = "percent mito genes",
        col = colb, las = 2, cex.axis = 0.7)

plot(qc[,"ribo_pct"], col = colb[batch], xlab = "cell index", 
     main = "% ribo (Rpl*) genes")
legend("topleft", legend = levels(batch), fill = colb, cex = 0.8)
boxplot(qc[,"ribo_pct"] ~ colData(se)$batch, main = "percent ribo genes",
        col = colb, las = 2, cex.axis = 0.7)
dev.off()

print(paste("Percent total variance captured in first 10 expression PCs is",
  round(cumsum(pca$sdev ^ 2 / sum(pca$sdev ^ 2))[10] * 100, digits = 2)
))

pdf(file = file.path(vizdir, pasteu(exptstr, "QCpca_prefilt.pdf")),
    height = 11, width = 8.5)
screeplot(pca, type = "lines", npcs = 50, main = "Expression PCA screeplot")
screeplot(qcpca, type = "lines", main = "QC-PCA screeplot")

fig_data <- data.frame(pca$x[, 1:2], qc, 
                       QPC1 = qcpca$x[, 1], QPC2 = qcpca$x[, 2],
                       batch = batch)

ggplot(fig_data, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(alpha = 0.3) + scale_color_manual(values = colb) +
  ggtitle("Expression PCA with simple filtering and scaling")
ggplot(fig_data, aes(x = QPC1, y = QPC2, color = batch)) +
  geom_point(alpha = 0.3) + scale_color_manual(values = colb) +
  ggtitle("QC metric PCA")

for (feature in
     c("log10_total_counts", 
       "ribo_pct",
       "mito_pct",
       "log10_total_features")) {
  pcaggplot(fig_data, "Expression PCA", feature)
  pcaggplot(fig_data, "QC metric PCA", feature)
}
dev.off()
 
  save(se, sce, qc, pca, qcpca,
       file = file.path(outdir, pasteu(exptstr, "prefilt.Rda")))
} else {
  load(file.path(outdir, pasteu(exptstr, "prefilt.Rda")))
}

# ---- Filtering ----
# store QC metrics with SummarizedExperiment object
colData(se) <- cbind(colData(se), qc)
se <- se[rowSums(assay(se)) > 0, ] # select expressed genes only

# select common genes
num_reads <- quantile(assay(se)[assay(se) > 0])[4]
num_cells <- 0.25 * ncol(se)
is_common <- rowSums(assay(se) >= num_reads) >= num_cells
table(is_common)

hk <- as.character(unlist(read.table(opt$hkfile)))
hk <- intersect(hk, rowData(se)$Symbol)

pdf(file = file.path(vizdir, pasteu(exptstr, "msf_prefilt.pdf")),
    title = "metric_sample_filtering")
mfilt <- metric_sample_filter(
  assay(se),
  nreads = colData(sce)$total_counts,
  gene_filter = is_common,
  pos_controls = hk,
  hard_nreads = 2000,
  zcut = 3,
  mixture = FALSE,
  plot = TRUE
)

msfpcapairsplot <- function(pca, mfilt, metric, title) {
  plot(
    pca$x,
    pch = 19,
    col = pal[as.numeric(mfilt[[metric]]) + 1],
    main = paste(toupper(deparse(substitute(pca))), "Filtered on", title)
  )
}

msfpcapairsplot(pca, mfilt, "filtered_breadth", "transcriptome 'breadth'")
msfpcapairsplot(pca, mfilt, "filtered_fnr", "FNR AUC")
msfpcapairsplot(pca, mfilt, "filtered_nreads", "nreads")
msfpcapairsplot(qcpca, mfilt, "filtered_breadth", "transcriptome 'breadth'")
msfpcapairsplot(qcpca, mfilt, "filtered_fnr", "FNR AUC")
msfpcapairsplot(qcpca, mfilt, "filtered_nreads", "nreads")

table(mfilt$filtered_nreads, mfilt$filtered_fnr)
filter_cell <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)

plot(qcpca$x,
     pch = 19,
     col = pal[as.numeric(filter_cell) + 1],
     main = "PCA Filtered")
dev.off()

# ---- Final gene filtering: highly expressed in at least N cells ----
num_reads <- quantile(assay(se)[assay(se) > 0])[4]
num_cells <- 3
is_quality <- rowSums(assay(se) >= num_reads) >= num_cells
table(is_quality)
se_filtered <- se[is_quality, filter_cell]
dim(se_filtered)
rownames(se_filtered) <- rowData(se_filtered)$Symbol
qc <- qc[colnames(se_filtered), ]
batch <- colData(se_filtered)$batch
expt <- colData(se_filtered)$expt

# Check for batch effects prior to normalization
pdf(file = file.path(vizdir, pasteu(exptstr, "batchchk_prefilt.pdf")))
plot(pca$x, pch = 19, col = colb[colData(se)$batch], 
     main = "PCA Color-coded by batch")
legend("topleft", legend = levels(batch), fill = colb, cex = 0.6)
plot(qcpca$x, pch = 19, col = colb[colData(se)$batch],
     main = "QPCA Color-coded by batch")
legend("bottomleft", legend = levels(batch), fill = colb, cex = 0.6)

boxplot(pca$x[, 1] ~ colData(se)$batch, col = colb, cex.axis = 0.75, las = 2,
        main = "First Principal Component")
boxplot(pca$x[, 2] ~ colData(se)$batch, col = colb, cex.axis = 0.75, las = 2,
        main = "Second Principal Component")

pca <- fastpca(log2(assay(se_filtered) + 1))
tsne_data <- Rtsne(pca[, 1:10], pca = FALSE, max_iter = 5000)
plot(tsne_data$Y, pch = 19, cex = 0.4, 
     col = alpha(colb[colData(se_filtered)$batch], 0.5))
legend("topleft", legend = levels(batch), fill = colb, cex = 0.5)
dev.off()

# ---- Preparation for normalization ----
poscon <- intersect(as.character(unlist(read.table(opt$posctrlfile))),
                    rowData(se_filtered)$Symbol)
rowData(se_filtered)$poscon <- rowData(se_filtered)$Symbol %in% poscon

# select negative controls (housekeeping)
hk <- intersect(hk, rowData(se_filtered)$Symbol)
negcon_eval <- sample(hk, length(poscon) / 2)
negcon_ruv <- setdiff(hk, negcon_eval)

rowData(se_filtered)$negcon_eval <- 
  rowData(se_filtered)$Symbol %in% negcon_eval
rowData(se_filtered)$negcon_ruv <- rowData(se_filtered)$Symbol %in% negcon_ruv

# ---- Gene expression plots post filtering, but pre-normalization ----
counts_filtered <- assay(se_filtered)
logcounts_filtered <- log2(counts_filtered + 1) 

controlheatmaps <- function(controllist, se_filtered){
  pdf(file = file.path(vizdir, pasteu(exptstr, controllist, "heatmap.pdf")))
  plotHeatmap(logcounts_filtered[rowData(se_filtered)[[controllist]], ],
              sampleData = data.frame(expt = colData(se_filtered)$expt, 
                                      batch = colData(se_filtered)$batch),
              clusterLegend = list(expt = cole, batch = colb),
              main = paste(toupper(controllist), "after filtering"), 
              breaks = .99)
  dev.off()
}
controlheatmaps("poscon", se_filtered)
controlheatmaps("negcon_eval", se_filtered)
controlheatmaps("negcon_ruv", se_filtered)

# Correlation of QC metrics with expression PCs for filtered data
cors <- sapply(1:5, function(i) abs(cor(pca[, i], qc, method = "spearman")))
bars <- data.frame(
  AbsoluteCorrelation = as.vector(cors),
  QC = factor(rep(colnames(qc), 5), levels = colnames(qc)),
  Dimension = as.factor(rep(paste0("PC", 1:5), each = ncol(qc)))
)

pdf(file = file.path(vizdir, pasteu(exptstr, "cor_qc_exprPCA.pdf")))
  ggplot(bars, aes(Dimension, AbsoluteCorrelation, group = QC, fill = QC)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_fill_manual(values = bigPalette) + ylim(0, 1) +
  ggtitle("Correlation between QC and expression PCA")
dev.off()

# ---- Save output ----
save(qc, pca, se_filtered, batch, expt,
     file = file.path(outdir, pasteu(exptstr, "se_filtered.Rda")))
save(counts_filtered, logcounts_filtered,
     file = file.path(outdir, pasteu(exptstr, "counts_filtered.Rda")))