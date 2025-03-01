#' Principle component analysis of all features color-coded by ED50 value prediction
#'
#' @description plotpca function generates a PCA plot showing all input features and highlighting dose-dependent features with color-coded ED50 values.
#' @usage plotpca(DoseResponse_report, DoseStat, EDrange = c(0,200),...)
#' @param DoseResponse_report list The output object of fitdrc() including ED50 values. See \code{\link{calcdosestat}}
#' @param DoseStat list The output object of calcdosestat(). See \code{\link{fitdrc}}
#' @param EDrange numeric The range of minimum and maximum ED50 value to be shown in the PCA plots.
#' @param ... Further arguments to be passed to plotpca.
#' @details plotPCA takes the raw intensity values and calculates the mean intensity for each feature for each dose. After this is calculated
#' range scaling is applied to make all intensities between zero and one. With the normalized mean intensities PCA is performed and the first
#' two components are plotted. Each dot represents a feature with the color representing the ED50 value. Grey dots were not determined to be
#' following a dose reponse trend.
#' @return plotpca function does not return any objects but generates a PDF with principal component plots of non-significant features (grey) and dose-response features with color-coded ED50 values.
#' @author Ethan Stancliff Yao <estancliff@wustl.edu>
#' Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' @import stats ggplot2 magrittr data.table
#'
plotpca = function(DoseResponse_report, DoseStat, EDrange = c(0,200),...){
  # Initialization
  if(!is.null(DoseStat)) {
    intensities <- DoseStat$Normalized_Response # obtain normalized responses from DoseStat
  } else {stop("DoseStat is not found.")}
  if(!is.null(DoseResponse_report$EDvalue)) {
    ED50s <- data.table(index=DoseResponse_report$EDvalue[,1], value = DoseResponse_report$EDvalue[,2])
    ED50s <- ED50s[value >= EDrange[1] & value<=EDrange[2]]  
  } else {stop("ED50 values is not found.")}
  if(!is.numeric(EDrange)||length(EDrange)!=2) {stop("EDrange must be specified as a pair of numeric values.")}
  cat("Plotting PCA and mapping ED50...")

  #make new dataframe of normalized mean intensities at each dose
  df = data.frame(intensities[,-1L])
  mat = data.frame(index = intensities[,1L])
  rep_index = 0
  for(replicate in DoseStat$Dose_Replicates){
    mat = cbind(mat,rowMeans(df[,(rep_index+1):(rep_index + replicate)]))
    rep_index = rep_index + replicate
  }
   colnames(mat) <- c("index",names(DoseStat$Dose_Replicates))
  #perform PCA
  pc2 = fortify(prcomp(mat[,-1L],scale. = TRUE,rank.=2))[,c("PC1","PC2")] %>% as.data.frame
  #get corresponding ED50 values
  tmp = as.numeric(rep(-1,nrow(pc2)))
  tmp[unlist(ED50s[,1])] <- unlist(ED50s[,2]) 
 
  #bind ED50 values to PCA
  pc2 = cbind(pc2, ED50=tmp)
  #plot the result
  p = ggplot2::ggplot() + geom_point(data=pc2[pc2$ED50 < -.5,], aes(x=PC1,y=PC2,color=ED50),alpha=.05) +
    geom_point(data=pc2[pc2$ED50 >= -.5,],aes(x=PC1,y=PC2,color=ED50),alpha=.5) +
    scale_color_gradient2(low="cyan", mid="blue",high="red",limits=c(0,max(pc2$ED50))) +
    theme_classic()
    #theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p
  filename = paste(paste(DoseResponse_report$projectName,"pca_ed50",sep="_"),".pdf",sep="")
  ggsave(filename = filename, width=8, height=8, plot = p, dpi = 300)
  cat("PCA plot is generated under: \n",getwd())
}

