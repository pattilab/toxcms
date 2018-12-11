plotDrc = function(dr4pl_Fit_result, doseResponse_report, mz_tag = "mzmed", rt_tag = "rtmed", file = "doseResponseCurve.pdf",...){

  # load dr4pl fitted result and corresponding feature data
  Feature <- doseResponse_report$Feature
  
  #index <- seq.int(nrow(Feature)) # generate index values
  mz <- Feature %>% dplyr::select(contains(mz_tag)) # get the m/z  from Feature
  rt <- Feature %>% dplyr::select(contains(rt_tag)) # get the rt from Feature
  
  if(length(mz)==1){
    cat("Succefully found",colnames(mz)," in feature table.\n")
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
  
  new_index <- sort(mz[[1]],index.return =TRUE)$ix
  last <- tail(new_index,1)
  pdf(file=file, width = 8.5, height = 11)
  temp_grid <- list()
  j=1
  
  for (i in new_index){
    
    dr_object <- dr4pl_Fit_result$dr4pl_object[[i]]
    if(is.null(dr_object)) {
      cat("Feature #", i,"has no fitting data. Continue plotting...\n")
      next
      }
    
    plot <- plotDr4pl(dr_object,dose_transform = TRUE, indices.outlier = TRUE)
    plot <- plot + labs(title=paste("index:", i, "  mz:", mz[i], "  rt:", rt[i])) + theme(plot.title=element_text(size=10, hjust = 0.5, face="bold", color="black", lineheight=1))
    temp_grid[[j]] <- plot
    j=j+1
    
      if(j==7 || i==last){
      do.call(grid.arrange,c(temp_grid,nrow=3,ncol=2,newpage=TRUE))
      temp_grid <- list()
      j=1
      }
  }
  
  cat("Plotting finished.\n")
  dev.off()
 
}


