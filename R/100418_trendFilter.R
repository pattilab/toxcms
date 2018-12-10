trendCalc = function(DoseStat, pval_cutoff){
  
  Dose_Levels <- DoseStat$Dose_Levels
  pval_indicator <- DoseStat$pvalue[,2:ncol(DoseStat$pvalue)] %>% apply(2,function(x){ifelse(x<=pval_cutoff,1,0)})
  pval_sum <- rowSums(pval_indicator)
  relchange_indicator <- DoseStat$relChange[,2:ncol(DoseStat$relChange)] %>% apply(2,function(x){ifelse(x<0,-1,1)}) #problematic when relchange = 0
  test_val  = rowSums(pval_indicator * relchange_indicator) # place ABSOLUTE in # Condition 2
  result = data.table(index = DoseStat$pvalue$index, pval_indicator = pval_sum, trend_indicator = test_val)
  
  return(result)
}

trendFilter = function(DoseStat, pval_cutoff = 0.05, pval_thres = 1,anova_cutoff = 0.05,trend = "both", relChange_cutoff = NULL, export = FALSE){
  
  trendCalc_result = trendCalc(DoseStat, pval_cutoff)
  index_pval <- which(trendCalc_result$pval_indicator>=pval_thres)
  # if anova p value is assigned, filter anova p value
  if (is.numeric(anova_cutoff)){
    index_aov <- which(DoseStat$aov_pvalue<=anova_cutoff)
    index_pval <- intersect(index_pval,index_aov)
  }
  # trend filter
  if (trend == "increase"){
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator & trendCalc_result$trend_indicator >0)
  } else if (trend == "decrease") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator & trendCalc_result$trend_indicator <0)
  } else if (trend == "both") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) == trendCalc_result$pval_indicator)
  } else if (trend == "reverse") {
    index_trend <- which(abs(trendCalc_result$trend_indicator) != trendCalc_result$pval_indicator)
    index <- intersect(index_trend, index_pval)
    doseResponse_report <- list(Feature = DoseStat$Feature[index,], Normalized_intensities = DoseStat$Normalized_intensities[index,], statistics = cbind(DoseStat$stat1[index,],DoseStat$stat2[index,]), trendCalc_result = trendCalc_result[index,], SampleInfo = DoseStat$SampleInfo, Dose_Replicates = DoseStat$Dose_Replicates, pval_cutoff=pval_cutoff,pval_thres = pval_thres,relChange_cutoff=NULL, trend = trend,anova_cutoff= anova_cutoff)
    cat("There are ", length(index)," features remaining.")
    if(export){
      file_name <- paste(Sys.Date(),"pval",pval_cutoff,"pval_pass#",pval_thres,"anova_cutoff",anova_cutoff,"trend",trend,sep = "_")
      write.csv(doseResponse_report$Feature,paste(file_name,"features.csv",sep="_"))
      write.csv(doseResponse_report$statistics,paste(file_name,"statistics.csv",sep="_"))
    }
    return(doseResponse_report)
  } else stop("trend should be set to either 'increase', 'decrease', 'both' or 'reverse'. ")
  
  if (is.numeric(relChange_cutoff)) {
    Dose_Levels <- DoseStat$Dose_Levels
    test_val <- trendCalc_result$trend_indicator
    pval_indicator_reverse <- DoseStat$pvalue[,2:ncol(DoseStat$pvalue)] %>% apply(2,function(x){ifelse(x<=pval_cutoff,0,1)})
    relchange_indicator_reverse <- DoseStat$relChange[,2:ncol(DoseStat$relChange)] %>% apply(2,function(x){ifelse(abs(x)<=relChange_cutoff,0,ifelse(x>0,1,-1))})
    test_val_reverse=pval_indicator_reverse * relchange_indicator_reverse
    test_val_reverse_mul=test_val_reverse %>% apply(2,function(x){x*test_val})
    test_val_reverse_mul_indicator=test_val_reverse_mul %>% apply(1,function(x)any(x<0))
    index_relChange_cutoff = which(test_val_reverse_mul_indicator == FALSE)
    index <- intersect(index_relChange_cutoff, intersect(index_trend, index_pval))
    doseResponse_report <- list(Feature = DoseStat$Feature[index,], Normalized_intensities = DoseStat$Normalized_intensities[index,], statistics = cbind(DoseStat$stat1[index,],DoseStat$stat2[index,]), trendCalc_result = trendCalc_result[index,], SampleInfo = DoseStat$SampleInfo, Dose_Replicates = DoseStat$Dose_Replicates, pval_cutoff=pval_cutoff,pval_thres = pval_thres,relChange_cutoff=relChange_cutoff, trend = trend,anova_cutoff= anova_cutoff)
    cat("There are ", length(index)," features remaining.")
  } else {
    cat("relChange_cutoff is not applied.\n")
    index <- intersect(index_trend, index_pval)
    doseResponse_report <- list(Feature = DoseStat$Feature[index,], Normalized_intensities = DoseStat$Normalized_intensities[index,], statistics = cbind(DoseStat$stat1[index,],DoseStat$stat2[index,]), trendCalc_result = trendCalc_result[index,], SampleInfo = DoseStat$SampleInfo, Dose_Replicates = DoseStat$Dose_Replicates, pval_cutoff=pval_cutoff,pval_thres = pval_thres,relChange_cutoff=NULL, trend = trend,anova_cutoff= anova_cutoff)
    cat("There are ", length(index)," features remaining.")
  }
  
  if(export){
    file_name <- paste(Sys.Date(),"pval",pval_cutoff,"pval_pass#",pval_thres,"relChange",relChange_cutoff,"anova_cutoff",anova_cutoff,"trend",trend,sep = "_")
    write.csv(doseResponse_report$Feature,paste(file_name,"features.csv",sep="_"))
    write.csv(doseResponse_report$statistics,paste(file_name,"statistics.csv",sep="_"))
  }
  
  return(doseResponse_report)
  
}

