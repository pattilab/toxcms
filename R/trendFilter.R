trendFilter = function(DoseStat, pval_cutoff = 0.05, pval_thres = 1,anova_cutoff = 0.05,trend =c("increase","decrease","mono","reverse","all"),
                       relChange_cutoff = NULL, export = FALSE){

  trend <- match.arg(trend)
  # trend calculations
  window <- ncol(DoseStat$statistic_multicomp)/2-1 # col size of pvalue and relchange
  pvalue <- DoseStat$statistic_multicomp[,3:(window+2)]
  relChange <- DoseStat$statistic_multicomp[,(window+3):(2*window+2)]
  pval_indicator <- pvalue %>% apply(2,function(x){ifelse(x<=pval_cutoff,1,0)})
  pval_sum <- rowSums(pval_indicator)
  relchange_indicator <- relChange %>% apply(2,function(x){ifelse(x<0,-1,1)}) #problematic when relchange = 0
  test_val  = rowSums(pval_indicator * relchange_indicator) # place ABSOLUTE in # Condition 2
  trendCalc_result <- data.table(index = DoseStat$statistic_multicomp$index, pval_indicator = pval_sum, trend_indicator = test_val)

  index_pval <- which(trendCalc_result$pval_indicator>=pval_thres)
  # if anova p value is assigned, filter anova p value
  if (is.numeric(anova_cutoff)){
    index_aov <- which(DoseStat$statistic_multicomp$aov_pvalue<=anova_cutoff)
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
    index <- intersect(index_trend, index_pval)
    cat("There are ", length(index)," features remaining.")
  }
  
  doseResponse_report <- list(Feature = DoseStat$Feature[index,], Normalized_response = DoseStat$Normalized_Response[index,], 
                              statistic_basic = DoseStat$statistic_basic[index,], statisic_multicomp = DoseStat$statistic_multicomp[index,], trendCalc_result = trendCalc_result[index,], 
                              SampleInfo = DoseStat$SampleInfo, Dose_Replicates = DoseStat$Dose_Replicates, pval_cutoff=pval_cutoff, pval_thres = pval_thres, relChange_cutoff=relChange_cutoff, 
                              trend = trend, anova_cutoff= anova_cutoff)

  if(export){
    file_name <- paste(Sys.Date(),"pval",pval_cutoff,"pval_pass#",pval_thres,"relChange",relChange_cutoff,"anova_cutoff",anova_cutoff,"trend",trend,sep = "_")
    write.csv(doseResponse_report$Feature,paste(file_name,"features.csv",sep="_"))
    write.csv(doseResponse_report$statistic_basic,paste(file_name,"statistic_basic.csv",sep="_"))
    write.csv(doseResponse_report$statistic_multicomp,paste(file_name,"statistic_multicomp.csv",sep="_"))
   }
  return(doseResponse_report)
}

