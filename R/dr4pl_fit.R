dr4pl_Fit <- function(doseResponse_report, ED=0.5, Dose_values,export = TRUE){
  
  if(!is.numeric(Dose_values)){
    stop("Dose input should be numeric.")
  }
  Dose_Replicates <- doseResponse_report$Dose_Replicates
  
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
  
  doses <- rep(Dose_values,times = Dose_Replicates)
  num_feature <- nrow(doseResponse_report$Feature)
  
  # initiation
  dr4pl_objects <- list(rep(NaN,num_feature))
  ED_values <- rep(NaN,num_feature)
  
  # looping dr4pl and ED calculation
  for (i in 1: num_feature) {
    vals <- as.numeric(doseResponse_report$Normalized_intensities[i,])
    tryCatch({
      dr4pl_objects[[i]] <- dr4pl(vals~doses)
      theta <- dr4pl_objects[[i]]$parameters
      ED_values[i] <- ObservedED (ED,theta)
    },error = function(e){
      cat("error occurs when fitting feature #",i,"\n")
    })
  }
  dr4pl_fit_result <- list(dr4pl_objects=dr4pl_objects,ED_values=ED_values)
  names(dr4pl_fit_result)[2]<-paste("ED",ED,sep="_") # rename ED_values to indicate ED used
  # combine ED values with original data
  if(export){
    file_name <- paste(Sys.Date(),"pval",doseResponse_report$pval_cutoff,"pval_pass#",doseResponse_report$pval_thres,"relChange",doseResponse_report$relChange_cutoff,"anova_cutoff",doseResponse_report$anova_cutoff,"trend",doseResponse_report$trend,sep = "_")
    a <- cbind(doseResponse_report$Feature,ED=dr4pl_fit_result[[2]])
    write.csv(a,paste(file_name,"features","with ED",ED,".csv",sep="_"))
  }
  return(dr4pl_fit_result)
}