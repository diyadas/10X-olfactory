# Helper functions for normalization of 10X data
# Authors: Diya Das and Rebecca Chance
# Last revised: Fri Mar  1 10:58:06 2019

library(ggplot2)
pal <- clusterExperiment::bigPalette
colb <- c("dodgerblue3", "darkviolet", "goldenrod", "chartreuse3", "cadetblue2", "magenta", "chocolate1", "forestgreen", "plum", "azure3", "chocolate4", "darkslategray", "cornsilk", "aquamarine3", "burlywood3", "darkblue", "gray45")
cole <- c("cornflowerblue", "darkgoldenrod", "darkorchid", "darkorange2", "deeppink3", "cadetblue1", "azure4", "darkslateblue", "darkolivegreen1", "antiquewhite2")
colRKC<- c("chartreuse3", "firebrick3", "chocolate4","slategray2","darkviolet", 
          "darkorange2","pink", "gold", "deepskyblue", "pink3", 
          "deeppink", "grey36", "royalblue3", "mediumorchid1","grey0", 
          "cyan2", "darkseagreen4")
t1 <- theme(plot.background = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),
                     axis.ticks = element_blank(), 
                     legend.background = element_blank(), 
                     axis.text.x = element_blank(), 
                     axis.text.y = element_blank(),
                     legend.key = element_rect(fill="white"), 
                     panel.border = element_rect(fill = NA,colour = "black"),
                     axis.line = element_blank(),
                     aspect.ratio = 1)


pasteu <- function(...) paste(..., sep = "_")
pasteu0 <- function(...) {
  tmp <- pasteu(...)
  stringi::stri_replace_last_regex(tmp, "__.|_.", ".")
}

fastpca <- function(expr, scale = FALSE) {
  k <- 50
  svd_raw <- rARPACK::svds(scale(t(expr), center = TRUE, scale = scale), 
                           k = k, nu = k, nv = 0)
  pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:k])
  return(pc_raw)
}

detectedgeneplot <- function(se, plot_title = ""){
  hist(colSums(assay(se)), breaks = 30,
    xlab = 'Number of UMI',
    main = paste(plot_title, "Number of UMI per sample"))
  hist(colSums(assay(se) > 0) / nrow(se), breaks = 30,
    xlab = "Proportion of detected genes",
    main = paste(plot_title, "Proportion of detected genes"))
  hist(colSums(assay(se) > 0), breaks = 30,
    xlab = "Number of detected genes",
    main = paste(plot_title, "Number of detected genes"))
  hist(colMeans(assay(se) == 0), breaks = 30,
    xlab = "Proportion of zeros",
    main = paste(plot_title, "Proportion of zeros"))
  
  boxplot(colSums(assay(se)) ~ colData(se)$batch,
    main = paste(plot_title, "Number of UMI per sample"),
    col = colb, las = 2)
  boxplot(colSums(assay(se) > 0) ~ colData(se)$batch,
    main = paste(plot_title, "Number of detected genes"),
    col = colb, las = 2)
  boxplot(colMeans(assay(se) == 0) ~ colData(se)$batch,
    main = paste(plot_title, "Proportion of zeros"),
    col = colb, las = 2)
}

pcaggplot <- function(fig_data, plot_type, feature){
  if (plot_type == "Expression PCA"){
    p <- ggplot(fig_data, aes(x = PC1, y = PC2, color = get(feature))) 
  } else if (plot_type == "QC metric PCA"){
    p <- ggplot(fig_data, aes(x = QPC1, y = QPC2, color = get(feature)))
  }
  return(p + geom_point(alpha = 0.3) + 
           scale_color_continuous(low = "blue", high = "yellow", 
                                  name = feature) + 
           ggtitle(plot_type))
}

rmblank <- function(input_filepath){
   staplr::remove_pages(1, input_filepath, output_filepath = "tmp.pdf")
   system(paste("mv tmp.pdf", input_filepath))
}

last_datfile <-	function(datdir, exptstr, method, norm, clusmethod){
   datfiles <<- list.files(path = datdir,
                           pattern = pasteu(exptstr, method, norm, clusmethod),
                           full.names = TRUE)
   datfile <- datfiles[length(datfiles)]
   print(paste("Loading this data file: ", datfile))
}

showPalette <- function (colPalette = bigPalette, which = NULL, cex = 1) {
  oldPar <- par(no.readonly = TRUE)
  wh <- which
  if (is.null(wh)) {
    wh <- seq_along(colPalette)
  }
  else {
    colPalette <- colPalette[wh]
  }
  n <- ceiling(sqrt(length(colPalette)))
  nblank <- n^2 - length(colPalette)
  xwide <- n
  yup <- n
  x1 <- rep(c(seq_len(xwide)) - 0.5, yup)
  x2 <- rep(c(seq_len(xwide)) + 0.5, yup)
  xtext <- rep(c(seq_len(xwide)), yup)
  ycolor1 <- rev(rep(seq(1, yup * 2, by = 2) - 0.5, each = xwide))
  ycolor2 <- rev(rep(seq(1, yup * 2, by = 2) + 0.5, each = xwide))
  ytext <- rev(rep(seq(2, yup * 2, by = 2) + 0.5, each = xwide))
  par(mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0.5, xwide + 0.5), ylim = c(0.5, (yup * 
                                                           2) + 0.5))
  rect(x1, ycolor1, x2, ycolor2, col = c(colPalette, rep("white", 
                                                         nblank)), border = FALSE)
  if (length(colPalette) > 100) {
    half <- ceiling(length(colPalette)/2)
    adj.text <- cbind(rep(0.5, half * 2), rep(c(0, 1), times = half))
    adj.text <- adj.text[seq_along(colPalette), ]
  } else {
    adj.text <- matrix(c(0.5, 0), nrow = length(colPalette), 
                       ncol = 2, byrow = TRUE)
  }
  for (i in seq_along(colPalette)) {
    text(xtext[i], ytext[i] - 1, colPalette[i], cex = cex, 
         adj = adj.text[i, ])
    if (length(colPalette) <= 100) 
      text(xtext[i], ytext[i] - 2, wh[i], cex = cex, adj = c(0.5, 
                                                             1))
  }
  par(oldPar)
}

recolorMassive <- function(ce){
  for (x in colnames(ce@clusterMatrix)){
    ce <- recolorClusters(ce, massivePalette[1:length(unique(ce@clusterMatrix[,x]))], whichCluster = x, "name")
  }
  return(ce)
}

recolorMassive_cl <- function(cl){
  for (x in colnames(cl@clusterMatrix)){
    cl <- recolorClusters(cl, massivePalette[1:length(unique(cl@clusterMatrix[,x]))], whichCluster = x, "name")
  }
  return(cl)
}

