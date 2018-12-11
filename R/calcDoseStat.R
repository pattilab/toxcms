calcDoseStat = function(Feature, Dose_Levels, padjust = NULL){
  # initiation
  Feature = data.table(Feature)
  num_feature <- nrow(Feature)
  Dose_Replicates <- rep(0,length(Dose_Levels))
  names(Dose_Replicates) <- Dose_Levels
  Response <- data.table(index = seq.int(num_feature))
  SampleNames <- character()
  
  # compute replicate number by name match
  for(DoseX in Dose_Levels) {
    ResponseX <- Feature %>% dplyr::select(contains(DoseX)) # get the raw response value from X
    if(length(ResponseX)==0){
      cat("Warning message: the following dose does not exist:", DoseX, "\n","Omit dose", DoseX, "and continue...","\n")
      Dose_Replicates <- Dose_Replicates[-which(Dose_Levels==DoseX)]
      Dose_Levels <- Dose_Levels[-which(Dose_Levels==DoseX)]
    } else {
      Dose_Replicates[DoseX] <- length(ResponseX)    # number of replicates in doseX
      Response = cbind(Response,ResponseX)           # data clean-up
    }
  }
  
  if(length(Dose_Levels)<3){
    cat("The number of doses is less than 3. NULL is returned.")
    return(NULL) 
  }
  
  SampleNames = colnames(Response[,2:length(Response)])
  SampleInfo = cbind(SampleNames, rep(Dose_Levels, times = Dose_Replicates)) #BUG FIXED
  # Data normalization -------------------------------------------------
  
  # Formula: (x-min)/(max-min)
  Maxint <- apply(Response[,2:length(Response)],1,max) # the exact number should be replaced by a sorting function in future 
  Minint <- apply(Response[,2:length(Response)],1,min) # same as above
  Delta <- Maxint - Minint # the range of intensity
  Normalized_Response <- apply(Response[,2:length(Response)],2,function(x) (x-Minint)/Delta) %>% data.table # run for each feature(row) 
  Normalized_Response[is.na(Normalized_Response)]<-0 # cy
  
  # Statistical Analysis Part I: Basic statistics ---------------------------------------------
  
  # Initiation
  stat1 <- data.table(index = Response$index)
  mean  <- matrix(nrow = num_feature, ncol = 0)
  std  <- matrix(nrow = num_feature, ncol = 0)
  cv <- matrix(nrow = num_feature, ncol = 0)
  
  # Statistical calculation
  
  for (DoseX in Dose_Levels) {
    
    ResponseX <- Response %>% dplyr::select(contains(DoseX))
    meanX <- apply(ResponseX, MARGIN = 1,mean) # calculate the mean in doseX
    stdX <- apply(ResponseX, MARGIN = 1,sd) # calculate the std in doseX
    cvX <- stdX / meanX # calculate the coefficient of variance
    
    # data clean-up
    mean <- cbind(mean,meanX)
    std <- cbind(std,stdX)
    cv <- cbind(cv, cvX)
  }
  
  # rename the columns
  colnames(mean) <- sapply(Dose_Levels, function(x) paste("mean",x,sep = ""))
  colnames(std) <- sapply(Dose_Levels, function(x) paste("std",x,sep = ""))
  colnames(cv) <- sapply(Dose_Levels, function(x) paste("cv",x, sep = ""))
  
  stat1 <- cbind(stat1, mean, std, cv) # FUNCTION OUTPUT
  
  # Statistical Analysis Part II: Calulate inter-group differences--------------------------------------------
  
  # if multigroup is FALSE, only perform adjacent group comparison
  # if multigroup is TRUE, adjust p value
  if (!is.null(padjust)){
    # initiation
    stat2 <- data.table(index =Response$index)
    relChange <- matrix(ncol = 0, nrow = num_feature) # Relative changes in between groups 
    colname <- character(length = 0)
    # compute relative change for all pairwise comparison
    for (i in 1:(length(Dose_Levels)-1)){
      end_1 <- sum(Dose_Replicates[1:i])
      start_1 <- end_1-as.numeric(Dose_Replicates[i])+1
      current <- Normalized_Response[,start_1:end_1]
      for(j in (i+1):length(Dose_Levels)){
        end_2 <- sum(Dose_Replicates[1:j])
        start_2 <- end_2-as.numeric(Dose_Replicates[j])+1
        higher <- Normalized_Response[,start_2:end_2]
        tempcol <- paste(Dose_Levels[i],"and",Dose_Levels[j])
        relChange <- cbind(relChange,rowMeans(higher)-rowMeans(current))
        colname <- c(colname,tempcol)
      }
    }
    ## rename relChange
    colnames(relChange) <- colname
    ## specify dose level for pvalue calculation
    dose_level <- factor(rep(Dose_Levels,time=Dose_Replicates),levels=Dose_Levels)
    ## calculate adjusted p value, replace NaN with Inf
    pvalue <- apply(Normalized_Response, MARGIN =1, function(x) {a=pairwise.t.test(as.numeric(x),dose_level,p.adj = padjust)$p.value; return(a[lower.tri(a,TRUE)])})
    pvalue <- t(pvalue)
    pvalue[is.nan(pvalue)] <- Inf
    ## rename pvalue
    colnames(pvalue) <- colname
    # store relChange and pvalue result
    relChange_result = cbind(stat2, relChange)
    pvalue_result = cbind(stat2,pvalue)
  } else {
    # initiation
    stat2 <- data.table(index =Response$index)
    relChange <- matrix(ncol = 0, nrow = num_feature) # Relative changes in between groups 
    pvalue <- matrix(ncol = 0, nrow = num_feature)
    end<-Dose_Replicates[1]
    colname <- character(length = 0)
    # for each dose-pair, e.g., 1_2, 2_3, etc.
    for (i in 1:(length(Dose_Levels)-1)){
      start <- end+1-Dose_Replicates[i]
      higher <- Normalized_Response[,(end+1):(end+Dose_Replicates[i+1])]
      current <- Normalized_Response[,start:end]
      relChange <- cbind(relChange,rowMeans(higher)-rowMeans(current))
      # calculate pvalue
      pvalX = apply(Normalized_Response, MARGIN =1, function(x) {t.test(x[start:end],x[(end+1):(end+Dose_Replicates[i+1])],var.equal = FALSE)$p.value})
      # data clean-up
      pvalue = cbind(pvalue, pvalX)
      end <- end+Dose_Replicates[i+1]
      tempcol <- paste(Dose_Levels[i],"and",Dose_Levels[i+1])
      colname <- c(colname,tempcol)
    }
    # replace pvalue equals NAN with Inf
    pvalue[is.nan(pvalue)] <- Inf # cy
    # rename the columns
    colnames(relChange) <- colname
    colnames(pvalue) <- colname
    # store relChange and pvalue result
    relChange_result = cbind(stat2, relChange)
    pvalue_result = cbind(stat2,pvalue)
  }
  ## specify dose level for pvalue calculation
  dose_level <- factor(rep(Dose_Levels,time=Dose_Replicates),levels=Dose_Levels)
  ## calculate anova p value
  aov_pvalue <- apply(Normalized_Response, MARGIN =1, function(x) {a=aov(as.numeric(x)~dose_level); return(as.numeric(unlist(summary(a))["Pr(>F)1"]))})
  # data final clean-up
  DoseStat = list(Feature = Feature, Dose_Levels = Dose_Levels, Normalized_intensities = Normalized_Response, stat = stat, relChange = relChange_result, pvalue = pvalue_result, Dose_Replicates = Dose_Replicates,aov_pvalue = aov_pvalue)
  return(DoseStat)
}
