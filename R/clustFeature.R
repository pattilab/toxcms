clustFeature <- function(DoseResponse_report, reference_index, sort.method = c("clust","range"), sort.thres = 20, 
                         dist.method = "euclidean", hclust.method = "average", mztag = "mzmed", rttag = "rtmed",
                         heatmap.on = FALSE, plot.all = FALSE,filename = "Heatmap.pdf",...){
  
  # choices are c("dist") c("clust","range") c("clust", "layer")
  # execute appropriate sorting method and return a list of indexes for each feature
  # plotting is finished here
  if (sort.method[1] == "dist"){
    cat("Extracting features based on distances...\n")
    if(sort.method[2]=="near"||sort.method[2] =="far"||sort.method[2]=="both"){
    index_to_extract <- sortByDist(DoseResponse_report, Index = reference_index, Sort.method = sort.method[2], Sort.thres = sort.thres, Dist.method = dist.method, Hclust.method = hclust.method, 
                                   mz_tag = mztag, rt_tag = rttag, Heatmap.on = heatmap.on, Plot.all = plot.all, Filename = filename )
    } else {stop("The 'sort.method' for dist is not specified. Should be 'near', 'far' or 'both'.")}
  } else if (sort.method[1] == "clust") {
    cat("Extracting features based on herarchical clustering...\n")
    if(sort.method[2]=="layer"||sort.method[2] =="range"){
    index_to_extract <- sortByClust(DoseResponse_report, Index = reference_index, Sort.method = sort.method[2], Sort.thres = sort.thres, Dist.method = dist.method, Hclust.method = hclust.method, 
                                    mz_tag = mztag, rt_tag = rttag, Heatmap.on= heatmap.on, Plot.all = plot.all, Filename = filename)
    } else { stop("The 'sort.method' for 'clust' is not specified. Should be 'layer' or 'range'. ") }
  } else { stop("The 'sort.method' do not exist.") }
  
  # extract the features based on the index
  FeatureClust = list()
  j=1
  for (i in reference_index){
    ind <- index_to_extract[[as.character(i)]]
    Feature_i <- cbind(New_index = ind, DoseResponse_report$Feature[ind,],DoseResponse_report$Normalized_Response[ind,-1L], DoseResponse_report$stat[ind,-1L], cbind(doseResponse_report$pvalue,doseResponse_report$relChange, doseResponse_report$trendCalc_result)[ind,-1L])
    FeatureClust[[j]] <- Feature_i
    j=j+1
  }
  # names are feature index
  names(FeatureClust) <- reference_index
  return(FeatureClust)
} 

sortByDist <- function(DoseResponse_report, Index, Sort.method=c("near","far","both"), Sort.thres = 10, Dist.method = "euclidean", Hclust.method = "average", 
                       mz_tag = mztag, rt_tag = rttag, Heatmap.on = FALSE,  Plot.all=FALSE, Filename = "Heatmap.pdf",...){
  
  Sort.method = match.arg(Sort.method)
  dataset <- DoseResponse_report$Normalized_response[,-1L]
  dist <- dataset %>% dist(method=Dist.method) %>% as.matrix
  
  # define the number of nearby features to extract and clustered
  if(Sort.thres<1 && Sort.thres>0){
    length_to_extract <- dist %>% nrow %>% `*`(Sort.thres) %>% ceiling
  } else if (Sort.thres>=1 && Sort.thres%%1 ==0 && Sort.thres<=nrow(dist)) {
    length_to_extract <- Sort.thres
  } else{ stop("'sort.thres' should be a fraction between 0 and 1 or an integer less than the number of features (",nrow(dist),").") }
  
  # store the extracted feature indexes
  
  Index_to_extract <- list()
  j=1
  # extract index
  if(Sort.method=="near"){
    extract_index <- 1:(length_to_extract+1)
  } else if(Sort.method=="far"){
    extract_index <- c(1,(ncol(dist)-length_to_extract+1):ncol(dist))
  } else if(Sort.method=="both"){
    extract_index <- c(1:(length_to_extract+1),(ncol(dist)-length_to_extract+1):ncol(dist))
  } else stop("Sort.method is not recognized.")
  
  for (i in Index){
  index_to_extract_i <- sort(dist[i,], index.return = TRUE)$ix[extract_index]
  # not efficient
  Index_to_extract[[j]] = index_to_extract_i
  j=j+1
  }
  
  names(Index_to_extract) <- Index
  
  # mztag and rttag search
  mz <- DoseResponse_report$Feature %>% dplyr::select(contains(mz_tag)) # get the m/z  from Feature
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
  
  
  # plotting the heatmap
  if(Heatmap.on==T){
    pdf(file = Filename, width = 8.5, height = 11)
    colLabel <- rep(names(DoseResponse_report$Dose_Replicates),times=DoseResponse_report$Dose_Replicates)
    if(Plot.all==T){
    heatmap.2(as.matrix(dataset), Colv = FALSE, margins = c(10,10), distfun = function(x) dist(x,method=Dist.method), hclustfun = function(x) hclust(x,method=Hclust.method), 
              main = paste("Hierarchical Clustering of ",nrow(dataset)," features",sep = ""), dendrogram = 'row',
              labCol = colLabel, trace="none",key=TRUE, keysize=1.3, key.title = "Scale")
    cat("Heatmap is generated under: \n",getwd())  
    } 
      for(i in Index){
        ind <- Index_to_extract[[as.character(i)]]
        mz_plot <- as.numeric(unlist(mz[ind,]))
        rt_plot <- as.numeric(unlist(rt[ind,]))
        labrow <- paste0(ind, "_", sprintf("%.4f", mz_plot),"_",rt_plot)
        heatmap.2(as.matrix(dataset[ind,]), Colv = FALSE, margins = c(10,10), distfun = function(x) dist(x,method=Dist.method), hclustfun = function(x) hclust(x,method=Hclust.method), 
                  main = paste("Feature #",i," mz=", mz[i], " rt=", rt[i], sep = ""), dendrogram = 'row', 
                  labCol = colLabel, labRow = labrow, trace="none",key=TRUE, keysize=1.2, key.title = "Scale")
      }
      cat("Heatmap is generated under: \n",getwd())  
      dev.off()
  }
  
  return(Index_to_extract)
}

sortByClust <- function(DoseResponse_report, Index, Sort.method = "layer", Sort.thres = 50, Dist.method="euclidean", Hclust.method = "average", 
                        mz_tag = mztag, rt_tag = rttag, Heatmap.on = FALSE, Plot.all=FALSE, Filename = "Heatmap.pdf",...){

  dataset <- DoseResponse_report$Normalized_Response[,-1L]
  clust <- dataset %>% dist(method=Dist.method) %>% hclust(method=Hclust.method)
  
  Index_to_extract = list()
  j=1
  if(Sort.method == "layer"){
    if(Sort.thres >=1 && Sort.thres%%1==0 && Sort.thres <= length(clust$order)){
      clust_Index <- cutree(clust, Sort.thres)
      for(i in Index){
        Index_to_extract_i <- which(clust_Index == clust_Index[i])
        Index_to_extract[[j]] <- Index_to_extract_i
        j=j+1
      }
      names(Index_to_extract) <- Index
    } else{ stop("When using 'clust-layer' sorting, the 'Sort.thres' should be a positive integer less than the number of features.") }
    
  } else if (Sort.method == "range"){
    if (Sort.thres <1 && Sort.thres>0) {
      length_to_extract <- clust$dist %>% as.matrix %>% nrow %>% `*`(Sort.thres) %>% ceiling
    } else if (Sort.thres >=1 && Sort.thres%%1==0 && Sort.thres <= length(clust$order)){ 
      length_to_extract <- Sort.thres
    } else { stop("When using 'clust-range' sorting, the 'Sort.thres' should be a fraction or a positive integer less than the number of features.")}
    tree <- cutree(clust,k=1:nrow(dataset))
    # for each index, find the target cluster that has at least `length_to_extract` features  
    for (i in Index) {
      numFeature_per_layer <- apply(tree,2,function(x) length(which(x==x[i])))
      # find the first layer from highest that contains no more than length_to_extract number of features
      target_layer <- min(which(numFeature_per_layer<=length_to_extract))
      clust_Index <- tree[,target_layer]
      Index_to_extract_i <- which(clust_Index == clust_Index[i])
      Index_to_extract[[j]] <- Index_to_extract_i
      j=j+1
    }
    names(Index_to_extract) <- Index
    
  } else { stop("The specified 'sort.method' of SortByClust() does not exist.") }
  
  # mztag and rttag search
  mz <- DoseResponse_report$Feature %>% dplyr::select(contains(mz_tag)) # get the m/z  from Feature
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
  
  if(Heatmap.on==T){
    pdf(file = Filename, width = 8.5, height = 11)
    colLabel <- rep(names(DoseResponse_report$Dose_Replicates),times=DoseResponse_report$Dose_Replicates)
    
    if(Plot.all==T){
      heatmap.2(as.matrix(dataset), Colv = FALSE, margins = c(10,10), distfun = function(x) dist(x,method=Dist.method), clustfun = function(x) hclust(x,method=Hclust.method), 
                main = paste("Hierarchical Clustering of ",nrow(dataset)," features",sep = ""), dendrogram = 'row',
                labCol = colLabel, trace="none",key=TRUE, keysize=1.3, key.title = "Scale")
    }
      for(i in Index){
        ind <- Index_to_extract[[as.character(i)]]
        mz_plot <- as.numeric(unlist(mz[ind,]))
        rt_plot <- as.numeric(unlist(rt[ind,]))
        labrow <- paste0(ind, "_", sprintf("%.4f", mz_plot),"_",rt_plot)
        heatmap.2(as.matrix(dataset[ind,]), Colv = FALSE, margins = c(10,10), distfun = function(x) dist(x,method=Dist.method), hclustfun = function(x) hclust(x,method=Hclust.method), 
                  main = paste("Feature #",i," mz=", mz[i], " rt=", rt[i], sep = ""), dendrogram = 'row', 
                  labCol = colLabel, labRow = labrow, trace="none",key=TRUE, keysize=1.3, key.title = "Scale")   
        }
    cat("Heatmap is generated under: \n",getwd(), "\n")  
    dev.off()
  }
  return(Index_to_extract)
}

# reference index refers to the row number in doseResponse_report

