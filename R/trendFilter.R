#' Trend filtering for metabolomic dose-response analysis (TOXcms step 2)
#'
#' @description The trendfilter() filters two pre-defined trends (dose-response models) from data processed by calcdosestat: monotonic ("increase","decrease","mono") trend and reverse trend.
#' The function returns the same list object as the input with only filtered trends retain. The resulting object is called doseResponse_report.
#' @usage trendfilter(DoseStat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05,trend =c("increase","decrease","mono","reverse","all"),
#' relChange_cutoff = NULL, export = FALSE)
#' @param DoseStat the output list object of calcdosesstat().
#' @param pval_cutoff a filtering threshold of significance cutoff of Welch t-test or post-hoc tests.
#' @param pval_thres a filtering threshold of numbers of significant dose pairs. With n dose levels, maximum value for adjacent comparison (multicomp="none") is n-1, maximum value for multi-group comparison
#' (multicomp="ttest"/"tukey"/"games-howell") is n(n-1)/2.
#' @param anova_cutoff a filtering threshold of significancce cutoff of ANOVA.
#' @param trend Type of trend. Options are "increase", "decrease", "mono" (both increase/decrease), "reverse" (one inflection point, "V" or "∧" shape), "all" (no filter).
#' @param relChange_cutoff a filtering threshold for monotonic trend ("increase"/"decrease"/"mono"). This filter removes the trends with non-significant relative change larger than the threshold and opposed to desired trends.
#' @param export a TRUE/FALSE value. When set to TRUE, a csv file including all results will be exported.
#' @return The function returns doseresponse_report, a list object including same elements as the input that are filtered based on trend and statistical thresholds.
#' @author Cong-Hui Yao <conghui.yao@wustl.edu>
#' Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' @import magrittr utils data.table
#' @export
trendfilter = function(DoseStat, pval_cutoff = 0.05, pval_thres = 1, anova_cutoff = 0.05,trend =c("increase","decrease","mono","reverse","all"),
                       relChange_cutoff = NULL, export = FALSE){

  trend <- match.arg(trend)
  # trend calculations
  pvalue <- DoseStat$pvalue[,c(-1L,-2L)]
  relChange <- DoseStat$relChange[,-1L]
  pval_indicator <- pvalue %>% apply(2,function(x){ifelse(x<=pval_cutoff,1,0)})
  pval_sum <- rowSums(pval_indicator)
  relchange_indicator <- relChange %>% apply(2,function(x){ifelse(x<0,-1,1)}) #problematic when relchange = 0
  test_val  = rowSums(pval_indicator * relchange_indicator) # place ABSOLUTE in # Condition 2
  trendCalc_result <- data.table(index = DoseStat$Feature$index, pval_indicator = pval_sum, trend_indicator = test_val)
  index_pval <- which(trendCalc_result$pval_indicator>=pval_thres)

  # if anova p value is assigned, filter anova p value
  if (is.numeric(anova_cutoff)){
    index_aov <- which(DoseStat$pvalue$aov_pvalue<=anova_cutoff)
    index_pval <- intersect(index_pval,index_aov)
  }
  # trend filter
  if (trend == "increase"){
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator & trendCalc_result$trend_indicator >0)
  } else if (trend == "decrease") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator & trendCalc_result$trend_indicator <0)
  } else if (trend == "mono") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator)
  } else if (trend == "reverse") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) != trendCalc_result$pval_indicator)
  } else if (trend == "all"){
    index_trend <- index_pval
  }

  if (is.numeric(relChange_cutoff) && !is.na(match(trend,c("increase","decrease","mono")))) {
    test_val <- trendCalc_result$trend_indicator
    pval_indicator_reverse <- pvalue %>% apply(2,function(x){ifelse(x<=pval_cutoff,0,1)})
    relchange_indicator_reverse <- relChange %>% apply(2,function(x){ifelse(abs(x)<=relChange_cutoff,0,ifelse(x>0,1,-1))})
    test_val_reverse=pval_indicator_reverse * relchange_indicator_reverse
    test_val_reverse_mul=test_val_reverse %>% apply(2,function(x){x*test_val})
    test_val_reverse_mul_indicator=test_val_reverse_mul %>% apply(1,function(x)any(x<0))
    index_relChange_cutoff = which(test_val_reverse_mul_indicator == FALSE)
    index <- intersect(index_relChange_cutoff, intersect(index_trend, index_pval))
    cat("There are ", length(index)," features remaining.")
  } else {
    cat("relChange_cutoff is not applied.\n")
    relChange_cutoff <- NULL
    index <- intersect(index_trend, index_pval)
    cat("There are ", length(index)," features remaining.")
  }

  parameters <- list(pval_cutoff=pval_cutoff, anova_cutoff= anova_cutoff, pval_thres = pval_thres, relChange_cutoff=relChange_cutoff, trend = trend)

  DoseResponse_report <- list(Feature = DoseStat$Feature[index,], Normalized_Response = DoseStat$Normalized_Response[index,],
                              stat = DoseStat$stat[index,], pvalue = DoseStat$pvalue[index,], relChange = DoseStat$relChange[index,],
                              trendCalc_result = trendCalc_result[index,], Dose_Levels = DoseStat$Dose_Levels, Dose_Replicates = DoseStat$Dose_Replicates,
                              SampleInfo = DoseStat$SampleInfo, projectName = DoseStat$projectName, parameters = parameters)

  if(export){
    file_name <- paste(DoseResponse_report$projectName, "pval",pval_cutoff,"trend",trend,"pval_pass#",pval_thres,"relChange",relChange_cutoff,"anova_cutoff",anova_cutoff,Sys.time(),sep = "_")
    write.csv(DoseResponse_report$Feature,paste(file_name,"features.csv",sep="_"))
    write.csv(DoseResponse_report$stat,paste(file_name,"statistic_basic.csv",sep="_"))
    write.csv(cbind(DoseResponse_report$pvalue,DoseResponse_report$relChange, DoseResponse_report$trendCalc_result),paste(file_name,"statistic_comparison.csv",sep="_"))
   }

  return(DoseResponse_report)
}

