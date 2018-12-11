calcDoseStat = function(Feature, Dose_Levels, Feature, Dose_Levels, 
                         multicomp = c("none","ttest","tukey","games-howell"), 
                         p.adjust.method=c("none","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),...){

  # Initialization ----------------------------------------------------------
  Feature = data.table(Feature)
  num_feature <- nrow(Feature)
  Response <- data.table(index = seq.int(num_feature))
  Dose_Replicates <- rep(0,length(Dose_Levels))
  names(Dose_Replicates) <- Dose_Levels
  SampleNames <- character()
  
  # match the assigned arguments(ensure the input is one of the choices) 
  norm.method <- match.arg(norm.method)
  multicomp <- match.arg(multicomp)
  p.adjust.method <- match.arg(p.adjust.method)

  # Statistical Analysis Part I: Basic statistics ---------------------------------------------
  stat1 <- data.table(index = Response$index)
  mean  <- matrix(nrow = num_feature, ncol = 0)
  std  <- matrix(nrow = num_feature, ncol = 0)
  cv <- matrix(nrow = num_feature, ncol = 0)
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
    meanX <- apply(ResponseX, MARGIN = 1,mean) # calculate the mean in doseX
    stdX <- apply(ResponseX, MARGIN = 1,sd) # calculate the std in doseX
    cvX <- stdX / meanX # calculate the coefficient of variance
    mean <- cbind(mean,meanX)
    std <- cbind(std,stdX)
    cv <- cbind(cv, cvX)
  }
  if(length(Dose_Levels)<3){
    cat("The number of doses is less than 3. NULL is returned.")
    return(NULL) 
  }
  colnames(mean) <- sapply(Dose_Levels, function(x) paste("mean",x,sep = ""))
  colnames(std) <- sapply(Dose_Levels, function(x) paste("std",x,sep = ""))
  colnames(cv) <- sapply(Dose_Levels, function(x) paste("cv",x, sep = ""))
  stat1 <- cbind(stat1, mean, std, cv) # FUNCTION OUTPUT
  SampleNames = colnames(Response[,2:length(Response)])
  SampleInfo = cbind(SampleNames, rep(Dose_Levels, times = Dose_Replicates)) #BUG FIXED
  
  # Data normalization -------------------------------------------------
  Normalized_Response <- normalization(Response[,2:length(Response)], Norm.method = norm.method) # "range' scaling by default
  Normalized_Response <- cbind(index=Response$index,Normalized_Response)
  Normalized_Response[is.na(Normalized_Response)]<-0 
  
  # Statistical Analysis Part II: Calulate inter-group differences--------------------------------------------
  stat2 <- data.table(index = Response$index)
  dose_levels <- factor(rep(Dose_Levels,times=Dose_Replicates),levels=Dose_Levels)
  dose_pair <- combn(unique(as.character(dose_levels)),2, paste,collapse="")

  # anova calculations
  aov_pvalue <- apply(Normalized_Response[,2:ncol(Normalized_Response)], MARGIN =1, function(x) {a=aov(as.numeric(x)~dose_levels); return(as.numeric(unlist(summary(a))["Pr(>F)1"]))})
      
  # multigroup test
   
  if(multicomp == "none" | multicomp =="ttest"){
  pairwise_pvalue <- apply(Normalized_Response[,2:ncol(Normalized_Response)], MARGIN =1, function(x) {
                  p <- pairwise.t.test(x,dose_levels, pool.sd=FALSE, p.adjust.method=p.adjust.method)$p.value
                  p <- p[lower.tri(p,diag=TRUE)]})
  } else{
    ind_tmp <- ifelse(multicomp=="tukey",1,2) # the first output table(1) is tukey, second(2) is games-howell
    pairwise_pvalue <- apply(Normalized_Response[,2:ncol(Normalized_Response)], MARGIN = 1, function(x){
                  p <- posthocTGH(x,dose_levels, method=multicomp, conf.level=0.95)$output[[ind_tmp]]$p}) # extract p value
  }
  pairwise_pvalue <- t(pairwise_pvalue)
  colnames(pairwise_pvalue) <- sapply(dose_pair,function(x) paste("p_",x,sep=""))
  pairwise_pvalue[is.nan(pairwise_pvalue)] <- 1 # maximum pvalue is 1.
  # calculate pairwise relative changes
  norm <- Normalized_Response[,2:ncol(Normalized_Response)]
  for (i in 1:(length(Dose_Levels)-1)){
      end_1 <- sum(Dose_Replicates[1:i])
      start_1 <- end_1-as.numeric(Dose_Replicates[i])+1
      current <- norm[,start_1:end_1]
      for(j in (i+1):length(Dose_Levels)){
        end_2 <- sum(Dose_Replicates[1:j])
        start_2 <- end_2-as.numeric(Dose_Replicates[j])+1
        higher <- norm[,start_2:end_2]
        relChange <- cbind(relChange,rowMeans(higher)-rowMeans(current))
      }
  }
  if(multicomp=="none"){
  dose_len <- length(Dose_Levels)
  seq_ind <- as.numeric()
  j <- 1
  for(i in 1:(dose_len-1)){
    ind <- j
    j <- j + dose_len-i
    seq_ind <- c(seq_ind,ind)
  }
  pairwise_pvalue<- pairwise_pvalue[,seq_ind]
  relChange <- relChange[,seq_ind]
  dose_pair <- dose_pair[seq_ind]
  }
  ## rename relChange
  colnames(relChange) <- sapply(dose_pair,function(x) paste("rel_",x,sep=""))
  # if adjacent pairwise applies
  stat2<- cbind(stat2,aov_pvalue,pairwise_pvalue,relChange)
  DoseStat = list(Feature = Feature, Normalized_intensities = Normalized_Response, statistical_basic = stat1, 
                  statistic_multicomp = stat2, Dose_Levels = Dose_Levels, Dose_Replicates = Dose_Replicates, 
                  SampleInfo = SampleInfo,norm.method="range",multicomp,p.adjust.method)
  return(DoseStat)
}

normalization <- function(data,Norm.method="range",...){
  data <- as.matrix(data)
  mean <- apply(data,MARGIN = 1,mean)
  std <- apply(data,MARGIN=1,sd)
  max <- apply(data,MARGIN=1,max)
  min <- apply(data,MARGIN=1,min)
  range <- max-min
  cv <- std/mean
  type = which(c("range","auto","pareto","vast","level","log10","log2","sqrt")==as.character(Norm.method))
  if(length(type)!=0){
     norm <- switch(type,
           apply(data,2,function(x) (x-min)/range), # range focus on correlation and restricted within 0 to 1
           apply(data,2,function(x) (x-mean)/std), # auto focus on correlation and restricted around 0
           apply(data,2,function(x) (x-mean)/sqrt(std)), # pareto focus on correlation and discriminate large fold-changes
           apply(data,2,function(x) (x-mean)/(std*cv)), # vast focus on less varying data and discriminate large varying data compared to auto scaling
           apply(data,2,function(x) (x-mean)/mean), # level scaling focus on fold changes
           data %>% log10 %>% apply(.,2,function(x) x-rowMeans(.)), # log10 transformation and center to mean, focus on large fold changes
           data %>% log2 %>% apply(.,2,function(x) x-rowMeans(.)), # log2
           data %>% sqrt %>% apply(.,2,function(x) x-rowMeans(.)) # sqrt pesudo scaling.
           )
     cat("Notice:",paste("'",Norm.method,"'",sep=""),"method applied for normalization.\n")
  } else stop(Norm.method, " is not a normalization method See ?calcDoseStat.")
  return(norm)
}


