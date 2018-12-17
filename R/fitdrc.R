fitdrc = function(DoseResponse_report, Dose_values, ED=0.5, mz_tag = "mzmed", rt_tag = "rtmed", export = TRUE, plot = TRUE, file = "doseResponseCurve.pdf",...){

  if(is.na(match(DoseResponse_report$parameters$trend,c("increase","decrease","mono")))){
    cat("fitdrc function works with monotonous trend analysis, namely the ‘trend’ parameter should be 'increase', 'decrease' or 'mono'.\n")
    cat("Suggestion: Use plottrend() instead.\n")
    stop("Aborted")
  }
  # load dr4pl fitted result and corresponding feature data
  cat("Fitting data into a logistic dose-response model 'dr4pl'...\n")
  dr4pl_Fit_result <- dr4pl_Fit(DoseResponse_report = DoseResponse_report, Dose_values = Dose_values, ED = ED, export = export)
  
  if(plot==TRUE){
    cat("\nPlotting fitted curves...\n")
  Feature <- DoseResponse_report$Feature
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
  dev.off()
  cat("Plotting finished. A pdf file is generated under:\n", getwd())
}

  ED_value <- cbind(DoseResponse_report$Feature$index,matrix(dr4pl_Fit_result[[2]],ncol=1))
  DoseResponse_report[["EDvalue"]] <- ED_value
  DoseResponse_report$parameters[["ED"]] <- ED
  return(DoseResponse_report)
}