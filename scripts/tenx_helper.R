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
pasteu0 <- function(...){
  tmp <- pasteu(...)
  stringi::stri_replace_last_fixed(tmp, "_.", ".")
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
                           pattern = pasteu(exptstr, method, norm, clusmethod), full.names = TRUE)
   datfile <- datfiles[length(datfiles)]
   print(paste("Loading this data file: ", datfile))
}