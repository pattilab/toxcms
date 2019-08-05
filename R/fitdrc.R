#'  Monotonic trend fitting and ED50 estimation for metabolomic dose-reponse analysis (TOXcms step 3)
#'
#' @description The fitdrc() fits the monotonic metabolic trends to a 4-parameter logistic model and visualize the fitting curve. Based on the fitted model,
#'     the function calculates the ED50 value for metabolic trend. The function returns the same list object as the input with additional table of ED50 values.
#'     The function also exports a PDF file including all plots. The output object is called doseResponse_report.
#' @usage fitdrc(DoseResponse_report, Dose_values, ED=0.5, mz_tag = "mzmed", rt_tag = "rtmed", export = TRUE, plot = TRUE,...)
#' @param DoseResponse_report list The output list object from trendfilter().
#' @param Dose_values numeric A vector of dose values.
#' @param ED numeric A value between (0,1), indicating the effective dose.
#' @param mz_tag character Name of the m/z column in the feature table.
#' @param rt_tag character Name of the retention time column in the feature table.
#' @param export TRUE/FALSE. When set to TRUE, a csv file including all results will be exported.
#' @param plot TRUE/FALSE. When set to TRUE, a PDF file including all plots will be exported.
#' @param ... Further arguments to be passed to fitdrc
#' @return The function returns doseresponse_report, a list object including same elements as the input and an addition of ED50 value.
#' @author Cong-Hui Yao <conghui.yao@wustl.edu>
#'     Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' @import magrittr dr4pl utils ggplot2 gridExtra grDevices
#' @export

fitdrc = function(DoseResponse_report, Dose_values, ED=0.5, mz_tag = "mzmed", rt_tag = "rtmed", export = TRUE, plot = TRUE,...){

  if(is.na(match(DoseResponse_report$parameters$trend,c("increase","decrease","mono")))){
    cat("fitdrc function works with monotonous trend analysis, namely the `trend` parameter should be 'increase', 'decrease' or 'mono'.\n")
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
  if(sapply(mz,is.numeric)) {
    mz <- round(mz,4)
  } else {
    mz %>% mutate_if(is.factor, as.character) %>% mz
  }
  if(sapply(rt,is.numeric)) {
    rt <- round(rt,4)
  } else {
    rt <- rt %>% mutate_if(is.factor, as.character) %>% data.table
  }
  new_index <- sort(dr4pl_Fit_result[[2]], index.return = TRUE)$ix
  last <- tail(new_index,1)
  pdf(file=paste(DoseResponse_report$projectName,"fitted_DoseResponseCurve", paste0("ED",round(ED*100)),".pdf",sep="_"), width = 8.5, height = 11)
  temp_grid <- list()
  j=1
  for (i in new_index){
    dr_object <- dr4pl_Fit_result$dr4pl_object[[i]]
    if(is.null(dr_object)) {
      cat("Feature #", i,"has no fitting data. Continue plotting...\n")
      next
      }
    plot <- plotDr4pl(dr_object,dose_transform = TRUE, indices.outlier = TRUE)
    plot <- plot + labs(title=paste0("index:", i, " ", mz_tag,":", mz[i], " ", rt_tag,":", rt[i,], paste0(" ED",round(ED*100),": ",round(dr4pl_Fit_result[[2]][i], 4), sep=""))) + theme(plot.title=element_text(size=10, hjust = 0.5, face="bold", color="black", lineheight=1))
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
