library(R.utils)
library(slingshot)

load("~/gitp/HBC-regen/output/clust/oeHBCregenWTKO/oeHBCregenWTKO_PCA.Rda")
col.palWT <- loadToEnv("~/gitp/HBC-regen/output/clust/oeHBCregenWT/oeHBCregenWT_PCA.Rda")[["col.pal"]]
clus.labelsWT <- loadToEnv("~/gitp/HBC-regen/output/clust/oeHBCregenWT/oeHBCregenWT_PCA.Rda")[["clus.labels"]]

XWT <- X[intersect(names(clus.labelsWT), rownames(X)),]
X2WT <- X2[intersect(names(clus.labelsWT), rownames(X2)),]
clus.labels2WT <- clus.labelsWT[rownames(X2WT)]
expt2WT <- expt2[rownames(X2WT)]
Eh="9"; Esus="3"; En="12";Ehrenew="1"
newss <- slingshot(X2WT[,1:6], clus.labels2WT, start.clus=Eh, end.clus=c(Esus, Ehrenew))

curves <- slingCurves(newss)
lineages <- slingLineages(newss)

library(rgl)
open3d()
plot3d(X[setdiff(rownames(XWT),rownames(X2WT)),1:3],col=col.pal[clus.labels[setdiff(rownames(XWT),rownames(X2WT))]],alpha=0.2, pch = 19, cex = 1, size=4, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
plot3d(X2WT[,1:3],col=col.palWT[clus.labels2WT],alpha=0.7, pch = 19, cex = 1, size=4, add=TRUE)
for (i in seq_along(curves)){
  plot3d(curves[[i]]$s[order(curves[[i]]$lambda),1:3],type="l",add=TRUE,lwd=2,col=col.palWT[which.max(tail(lineages[[i]],1)==levels(clus.labelsWT))])
}

X2KO <- X2[grep("SOX2CKO",expt2),]
clus.labels2KO <- clus.labels2[rownames(X2KO)]
expt2KO <- expt2[rownames(X2KO)]
newss_KO <- predict(newss, newdata=X2KO[,1:6])


for (i in 1:length(curves)){
  linedf <- data.frame(pseudotime = slingPseudotime(newss)[ ,i], lambda = curves[[i]]$lambda, w = curves[[i]]$w, clus.labels = clus.labels2WT, samples=rownames(X2WT), expt=expt2WT)
  #linedf$KO = linedf$expt %in% levels(linedf$expt)[7:12]
  linedf <- linedf[with(linedf, order(pseudotime)), ]
  
  medoids <- sapply(levels(linedf$clus.labels),function(clID){
    x.sub <- linedf$pseudotime[linedf$clus.labels == clID]
    col <- col.palWT[linedf$clus.labels][which.max(linedf$clus.labels==clID)]
    return(list(means = mean(x.sub, na.rm=TRUE), sdev= sd(x.sub, na.rm=TRUE), col=col))
  })
  means = unlist(medoids["means",])
  sdev = unlist(medoids["sdev",])
  col = unlist(medoids["col",])
  
  #svg(file.path(viz_dir, paste0(esh, "_Lineage", i, "_shrink", 1, "_", Sys.Date(),".svg")),width=8.5, height=2)
  par(mfrow=c(3, 1),mar=c(1, 1, 1, 1))
  plot(linedf$pseudotime,rep(0, length(linedf$pseudotime)),cex=3,axes=F, pch=16, xlab='', ylab='', col=alpha(col.palWT[linedf$clus.labels],0.5), xlim=c(min(linedf$pseudotime, na.rm=TRUE),max(linedf$pseudotime, na.rm=TRUE)),ylim=c(-0.1, 0.1)); abline(h=0, col=alpha("black", 0.5))
  points(x=means,y=rep(0.07, length(means)), col=col, pch=19)
  arrows(means-sdev, rep(0.07, length(means)), means+sdev, rep(0.07, length(means)), length=0.05, angle=90, code=3, col=col)
  #legend("left", levels(linedf$clus.labels), fill=col, cex=0.5, xpd=TRUE, inset=c(-0.12,0.1))
  
  
  plot(linedf$pseudotime,rep(0, length(linedf$pseudotime)),cex=3,axes=F, pch=16, xlab='', ylab='',col=alpha(cole[linedf$expt], 0.3),xlim=c(min(linedf$pseudotime, na.rm=TRUE),max(linedf$pseudotime, na.rm=TRUE)), ylim=c(-0.1, 0.1)); abline(h=0, col=alpha("black", 0.5))
  
  linedfKO <- data.frame(pseudotime = slingPseudotime(newss_KO)[ ,i], lambda = slingCurves(newss_KO)[[i]]$lambda, w = slingCurves(newss_KO)[[i]]$w, clus.labels = clus.labels2KO, samples=rownames(X2KO), expt=expt2KO)
  linedfKO <- linedfKO[with(linedfKO, order(pseudotime)), ]
  
  
  plot(linedfKO$pseudotime,rep(0, length(linedfKO$pseudotime)),cex=3,axes=F, pch=16, xlab='', ylab='',col=alpha(cole[linedfKO$expt], 0.3),xlim=c(min(linedf$pseudotime, na.rm=TRUE),max(linedf$pseudotime, na.rm=TRUE)), ylim=c(-0.1, 0.1)); abline(h=0, col=alpha("black", 0.5))
  
  #dev.off()
  
}
