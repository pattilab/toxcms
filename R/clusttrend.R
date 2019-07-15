#' Clustering trend for metabolomic dose-response analysis (TOXcms step 3)
#'
#' @description This function applies similarity matrices and hierarchical clustering to group metabolomic trends based on the similarity and dissimilarity in their dose-dependent responses.
#' @usage (DoseResponse_report, reference_index=NULL, sort.method = c("clust","layer"), sort.thres = 20,
#' dist.method = "euclidean", hclust.method = "average", mz_tag = "mzmed", rt_tag = "rtmed",
#' heatmap.on = FALSE, plot.all = FALSE,export=TRUE, filename = "Heatmap.pdf",...)
#' @param DoseResponse_report This is the output of trendfilter(), the end of toxcms step 2.
#' @param reference_index This is the indices of the reference trends. clusttrend() apply the algorithm for each of the reference index.
#' @param sort.method a two-element vector indicating the sorting methods. The first element can be either "dist" or "clust". If the first element is "dist", the second element can be "far","near" or "both".
#' If the first element is "clust", then the second-element can be “range" or "layer".
#' @param sort.thres a numeric number, indicating the sorting threshold for each group. Can be a fraction number smaller than 1 or a positibe integer larger than 1.
#' @param dist.method  Matrice used for similarity calculation. See also \link[stats]{dist}.
#' @param hclust.method Methods for hierarchical clustering. See also \link[stats]{hlust}.
#' @param mz_tag name of the m/z column in the feature table.
#' @param rt_tag name of the retention time column in the feature table.
#' @param heatmap.on a TRUE/FALSE value indicating whether heatmap shoudld be generated.
#' @param plot.all a TRUE/FALSE value indicating whether a heatmap of all trends should be generated.
#' @param export a TRUE/FALSE value indicating whetehr the trend clusters should be exported to a csv. file.
#' @param filename a character indicating the name of the exported csv file.
#' @param ... Further arguments to be passed to clusttrend.
#' @return clusttrend returns a list object, within wich each element is a trend cluster. The name of each element is the reference_index. Heatmap visualization of all trend clusters will be generated if heatmap.on=TRUE.
#' @author Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' Cong-Hui Yao <conghui.yao@wustl.edu>
#' @import ggplot2 magrittr grDevices data.table gridExtra dplyr
#' @importFrom stats hclust dist cutree
#' @importFrom gplots heatmap.2
#' @export
clusttrend <- function(DoseResponse_report, reference_index=NULL, sort.method = c("clust","layer"), sort.thres = 20,
                         dist.method = "euclidean", hclust.method = "average", mz_tag = "mzmed", rt_tag = "rtmed",
                         heatmap.on = FALSE, plot.all = FALSE,export=TRUE, filename = "Heatmap.pdf",...){

  # sort.method are c("dist","both"/"far"/"near") c("clust","range"/"layer")
  if(!is.null(reference_index)) {

    if (sort.method[1] == "dist"){
    cat("Sorting feature trends by distances...\n")
      if(sort.method[2]=="near"||sort.method[2] =="far"||sort.method[2]=="both"){
      index_to_extract <- sortByDist(DoseResponse_report, Index = reference_index, Sort.method = sort.method[2],
                                     Sort.thres = sort.thres, Dist.method = dist.method, Hclust.method = hclust.method)
      } else {stop("The 'sort.method' for dist is not specified. Should be 'near', 'far' or 'both'.")}
    } else if (sort.method[1] == "clust") {
    cat("Sorting feature trends by herarchical clustering...\n")
      if(sort.method[2]=="layer"||sort.method[2] =="range"){
      index_to_extract <- sortByClust(DoseResponse_report, Index = reference_index, Sort.method = sort.method[2],
                                      Sort.thres = sort.thres, Dist.method = dist.method, Hclust.method = hclust.method)
      } else { stop("The 'sort.method' for 'clust' is not specified. Should be 'layer' or 'range'. ") }
    } else { stop("The 'sort.method' do not exist.") }
  } else {
  cat("'reference index' is NULL. Sort all trends by hierarchical clustering...\n ")
    if(!(sort.method[1]=="clust"&&sort.method[2]=="layer")) {
    warning("When sorting all trends,'sort.method' are enforced to be 'clust-layer'. ")
    sort.method = c("clust","layer")
    }
  index_to_extract <- sortByClust(DoseResponse_report, Index = reference_index, Sort.method = sort.method[2],
                                  Sort.thres = sort.thres, Dist.method = dist.method, Hclust.method = hclust.method)
  }
# indexes of reference features for each cluster
ref_index <- as.numeric(names(index_to_extract))
FeatureClust = list() # sort features based on clusters
j=1
  for (i in ref_index){
  ind <- index_to_extract[[as.character(i)]]
  Feature_i <- cbind(clusterid=i,DoseResponse_report$Feature[ind,],DoseResponse_report$Normalized_Response[ind,-1L], DoseResponse_report$stat[ind,-1L], cbind(DoseResponse_report$pvalue,DoseResponse_report$relChange, DoseResponse_report$trendCalc_result)[ind,-1L])
  FeatureClust[[j]] <- Feature_i
  j=j+1
  }
# names are feature index
names(FeatureClust) <- names(index_to_extract)

# mz_tag and rt_tag search
mz <- DoseResponse_report$Feature %>% dplyr::select(contains(mz_tag)) # get the m/z from Feature
rt <- DoseResponse_report$Feature %>% dplyr::select(contains(rt_tag)) # get the rt from Feature

  if(length(mz)==1){
    cat("Succefully found",colnames(mz),"in feature table.\n")
  } else if(length(mz)==0) {
    stop("The specified '",mz_tag, "' does not match any columns in feature table.")
  } else if (length(mz)>1){
    stop("Multiple columns including '",mz_tag, "'. Please check the column names of feature table." )
  } else{
    stop("Error occurs when matching ", mz_tag)
  }
  if(length(rt)==1){
    cat("Succefully found",colnames(rt),"in feature table.\n")
  } else if(length(rt)==0){
    stop("The specified '",rt_tag, "'does not match any columns in feature table.")
  } else if (length(rt)>1){
    stop("Multiple columns including '",rt_tag, "'. Please check the column names of feature table." )
  } else{
    stop("Error occurs when matching ", rt_tag)
  }
mz <- round(mz,4)
rt <- round(rt)
dataset <- DoseResponse_report$Normalized_Response[,-1L]
  # plot the heatmap

  if(heatmap.on==T){
    pdf(file = filename, width = 8.5, height = 11)
    colLabel <- rep(names(DoseResponse_report$Dose_Replicates),times=DoseResponse_report$Dose_Replicates)
    if(plot.all==T){
    gplots::heatmap.2(as.matrix(dataset), Colv = FALSE, margins = c(10,10), distfun = function(x) stats::dist(x,method=dist.method), hclustfun = function(x) stats::hclust(x,method=hclust.method),
              main = paste("Hierarchical Clustering of ",nrow(dataset)," features",sep = ""), dendrogram = 'row',
              labCol = colLabel, trace="none", key=TRUE, keysize=1.3, key.title = "Scale")
    }
    for(i in ref_index){
    ind <- index_to_extract[[as.character(i)]]
    mz_plot <- as.numeric(unlist(mz[ind,]))
    rt_plot <- as.numeric(unlist(rt[ind,]))
    labrow <- paste0(ind, "_", sprintf("%.4f", mz_plot),"_",rt_plot)
    gplots::heatmap.2(as.matrix(dataset[ind,]), Colv = FALSE, margins = c(10,10), distfun = function(x) stats::dist(x,method=dist.method), hclustfun = function(x) stats::hclust(x,method=hclust.method),
              main = paste("Feature #",i," mz=", mz[i], " rt=", rt[i], sep = ""), dendrogram = 'row',
              labCol = colLabel, labRow = labrow, trace="none",key=TRUE, keysize=1.2, key.title = "Scale")
    }
  cat("Heatmap is generated under: \n",getwd())
  dev.off()
  }
return(FeatureClust)
}

sortByDist <- function(DoseResponse_report, Index, Sort.method=c("near","far","both"),
                       Sort.thres = 20, Dist.method = "euclidean", Hclust.method = "average",...){

Sort.method = match.arg(Sort.method)
dataset <- DoseResponse_report$Normalized_Response[,-1L]
dist <- dataset %>% stats::dist(method=Dist.method) %>% as.matrix
# calculate number of features to be exrtacted and clustered
  if(Sort.thres<1 && Sort.thres>0){
    length_to_extract <- dist %>% nrow %>% `*`(Sort.thres) %>% ceiling
  } else if (Sort.thres>=1 && Sort.thres%%1 ==0 && Sort.thres<=nrow(dist)) {
    length_to_extract <- Sort.thres
  } else{ stop("'sort.thres' must be a fraction between 0 and 1 or an integer less than the number of features (",nrow(dist),").") }

Index_to_extract <- list()
.j=1
  # extract indexes
  if(Sort.method=="near"){
    extract_index <- 1:(length_to_extract+1)
  } else if(Sort.method=="far"){
    extract_index <- c(1,(ncol(dist)-length_to_extract+1):ncol(dist))
  } else if(Sort.method=="both"){
    extract_index <- c(1:(length_to_extract+1),(ncol(dist)-length_to_extract+1):ncol(dist))
  } else stop("Sort.method is not recognized.")

  for (i in Index){
  index_to_extract_i <- sort(dist[i,], index.return = TRUE)$ix[extract_index]
  Index_to_extract[[.j]] = index_to_extract_i
  .j=.j+1
  }
names(Index_to_extract) <- as.character(Index)
return(Index_to_extract)
}

sortByClust <- function(DoseResponse_report, Index, Sort.method = c("layer","range"),
                        Sort.thres = 50, Dist.method="euclidean", Hclust.method = "average",...){
Sort.method = match.arg(Sort.method)
dataset <- DoseResponse_report$Normalized_Response[,-1L]
clust <- dataset %>% stats::dist(method=Dist.method) %>% stats::hclust(method=Hclust.method)
Index_to_extract = list()
.j=1

  if(is.null(Index)){
  Index = numeric()
    if(Sort.thres >=1 && Sort.thres%%1==0 && Sort.thres <= length(clust$order)){
    cat(nrow(dataset),"trends will be sorted into",Sort.thres,"clusters...\n")
    } else {
    Warning("'sort.thres' must be a positive integer no larger than number of features. Set to 50 by default and continute...")
    Sort.thres=50
    }
  clust_Index <- stats::cutree(clust, k=Sort.thres)
    for( .ind in unique(clust_Index)){
    Index_to_extract_i <- which(clust_Index == .ind)
      if(length(Index_to_extract_i)==1) next;
    Index_to_extract[[.j]] <- Index_to_extract_i
    Index[.j] <- Index_to_extract_i[1]
    .j=.j+1
    }
  names(Index_to_extract) <- as.character(Index)
  } else {
    if(Sort.method == "layer"){
      if(Sort.thres >=1 && Sort.thres%%1==0 && Sort.thres <= length(clust$order)){
      cat(nrow(dataset),"trends will be sorted into",Sort.thres,"clusters...\n")
      clust_Index <- stats::cutree(clust, k=Sort.thres)
        for(i in Index){
        Index_to_extract_i <- which(clust_Index == clust_Index[i])
        Index_to_extract[[.j]] <- Index_to_extract_i
        .j=.j+1
        }
      } else{ stop("When using 'clust-layer' method, 'sort.thres' must be a positive integer no larger than number of features.") }
    } else if (Sort.method == "range"){
      if (Sort.thres <1 && Sort.thres>0) {
      length_to_extract <- clust$dist %>% as.matrix %>% nrow %>% `*`(Sort.thres) %>% ceiling
      } else if (Sort.thres >=1 && Sort.thres%%1==0 && Sort.thres <= length(clust$order)){
      length_to_extract <- Sort.thres
      } else { stop("When using 'clust-range' method, the 'sort.thres' must be a fraction or a positive integer less than number of features.")}
    tree <- stats::cutree(clust,k=1:nrow(dataset))
    # for each index, find the target cluster that has at least `length_to_extract` features
      for (i in Index) {
      numFeature_per_layer <- apply(tree,2,function(x) length(which(x==x[i])))
      # find the first layer from highest that contains no more than length_to_extract number of features
      target_layer <- min(which(numFeature_per_layer<=length_to_extract))
      clust_Index <- tree[,target_layer]
      Index_to_extract_i <- which(clust_Index == clust_Index[i])
      Index_to_extract[[.j]] <- Index_to_extract_i
      .j=.j+1
      }
    } else { stop("The specified 'sort.method' does not exist.") }
  names(Index_to_extract) <- as.character(Index)
  }
return(Index_to_extract)
}
