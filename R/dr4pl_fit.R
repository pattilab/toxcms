#'  Fitting dose response four parameter logistic function to dose-dependent data points
#'
#' @description This function applies dr4pl function to fitting a dose-response four parameter logistic model to 
#' each dose-dependent metabolomic features. 
#' @usage  dr4pl_fit(DoseResponse_report, ED=0.5, Dose_values, export = TRUE)
#' @param DoseResponse_report list a complex list of dose-dependent features processed by trendfilter function.
#' @param ED numeric effective dose value from 0 to 1. 0.5 by default.
#' @param Dose_values dose values to be fitted.
#' @param export TRUE/FALSE whether to export the fitting results to a CSV file.
#' @import dr4pl

dr4pl_fit <- function(DoseResponse_report, ED=0.5, Dose_values, export = TRUE){

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
      dr4pl_objects[[i]] <- dr4pl(vals~doses, method.init="logistic")
      theta <- dr4pl_objects[[i]]$parameters
      #ED_values[i] <- exp(theta[2])
      ED_values[i] <- ObservedED(max(vals)*ED,theta)
    },error = function(e){
      cat("error occurs when fitting feature #",i,"\n")
    })
  }
  ED_values[is.nan(ED_values) | is.na(ED_values) | is.null(ED_values)] <- 0
  dr4pl_fit_result <- list(dr4pl_objects=dr4pl_objects,ED_values=ED_values)
  names(dr4pl_fit_result)[2]<-paste("ED",ED,sep="_") # rename ED_values to indicate ED used
  # combine ED values with original data
  if(export){
    param <- DoseResponse_report$parameters
    file_name <- paste(DoseResponse_report$projectName, "trend", param$trend, "pval",param$pval_cutoff,"pval_pass#",param$pval_thres,"relChange",param$relChange_cutoff, "anova_cutoff",param$anova_cutoff, paste0("ED",round(ED*100)), sep = "_")
    a <- cbind(DoseResponse_report$Feature,ED=dr4pl_fit_result[[2]])
    write.csv(x = a,file = paste(file_name,"features","with ED",ED,".csv",sep="_"))
    cat("\n Fitting result is exported under:\n",getwd())
  }
  return(dr4pl_fit_result)
}
