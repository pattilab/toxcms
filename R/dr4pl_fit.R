dr4pl_Fit <- function(DoseResponse_report, ED=0.5, Dose_values, export = TRUE){

  if(!is.numeric(Dose_values)){
    stop("Dose input should be numeric.")
  }
  Dose_Replicates <- DoseResponse_report$Dose_Replicates

  if(length(Dose_values)!=length(Dose_Replicates)){
    cat("Dose values",Dose_values,"do not match",names(Dose_Replicates),".\n")
    stop("Abort")
  }

  # print Dose values and corresponding Dose names
  cat("Dose values",Dose_values,"are assigned to",names(Dose_Replicates),".\n")

  # invalid ED input is forced to be 0.5
  if(ED<0 | ED>1){
    cat("Input ED is invalid.","\n ED_0.5 is returned.")
    ED <- 0.5
  }

  doses <- rep(Dose_values,times = Dose_Replicates)
  num_feature <- nrow(DoseResponse_report$Feature)

  # initiation
  dr4pl_objects <- list()
  ED_values <- rep(NaN,num_feature)

  # looping dr4pl and ED calculation
  for (i in 1: num_feature) {
    vals <- as.numeric(DoseResponse_report$Normalized_Response[i,-1L])
    tryCatch({
      dr4pl_objects[[i]] <- dr4pl(vals~doses)
      theta <- dr4pl_objects[[i]]$parameters
      ED_values[i] <- ObservedED(ED,theta)
    },error = function(e){
      cat("error occurs when fitting feature #",i,"\n")
    })
  }
  ED_values[is.nan(ED_values) | is.na(ED_values) | is.null(ED_values)] <- 0
  dr4pl_fit_result <- list(dr4pl_objects=dr4pl_objects,ED_values=ED_values)
  names(dr4pl_fit_result)[2]<-paste("ED",ED,sep="_") # rename ED_values to indicate ED used
  # combine ED values with original data
  if(export){
    file_name <- paste(DoseResponse_report$projectName, "trend",DoseResponse_report$parameters$trend,"pval",DoseResponse_report$parameters$pval_cutoff,"pvalthres#",DoseResponse_report$parameters$pval_thres,"relChange",DoseResponse_report$parameters$relChange_cutoff,"anova_cutoff",DoseResponse_report$parameters$anova_cutoff,sep = "_")
    a <- cbind(DoseResponse_report$Feature,ED=dr4pl_fit_result[[2]])
    write.csv(x = a,file = paste(file_name,"features","with ED",ED,".csv",sep="_"))
    cat("\n Fitting result is exported under:\n",getwd())
  }
  return(dr4pl_fit_result)
}

# function that calculate observed ED according to fitted parameters
  ObservedED <- function (ED, theta) {
    if(any(is.na(theta))) {
      stop("One of the parameter values is NA.")
    }
    if(theta[2]<=0) {
      stop("An IC50/ED50 estimate should always be positive.")
    }
    f <- as.numeric(theta[2]*((theta[4]-theta[1])/(ED-theta[1])-1)^(1/theta[3]))
    return(f)
  }
