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
