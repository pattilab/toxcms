#' Statistical analysis for metabolomic dose-response analysis (TOXcms step 1)
#'
#' @description The function calcdosestat() performs statistical analysis on metabolic responses at diffrent dose levels. Basic statistics of each dose level is calculated.
#' MS raw intensities of each feature are normalized by range scaling (x-min/max-min). ANOVA and multigroup testing are performed to obtain p-values and relative changes
#' among dose levels. The output is called dosestat.
#' @usage calcdosestat(Feature, Dose_Levels, multicomp = c("none","ttest","tukey","games-howell"),
#' p.adjust.method = c("none","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),...)
#' @param Feature data.table a feature table including basic information (m/z, retention time, etc) and MS raw intensities among all dose levels.
#' @param Dose_levels character Keywords indicating ordered dose levels in feature table. The dose levels should be ordered incrementally.
#' @param multicomp options for multigroup testing of significant difference. "none"- Welch t-test on adjacent doses; "ttest" - pairwise Welch t-test among all doses;
#' "tukey"/"games-howell" - post-hoc tests followed by ANOVA.
#' @param p.adjust.method p-value adjustment for false-positive reduction in multiple testing. See also \link[stats]{p.adjust}
#' @param ... Further arguments to be passed to calcdosestat.
#' @return calcdosestat returns a list object consisting all the statistical results described above.
#' @author Cong-Hui Yao <conghui.yao@wustl.edu>
#' Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' @import magrittr data.table dplyr
#' @importFrom stats aov pairwise.t.test sd
#' @importFrom userfriendlyscience posthocTGH
#' @export

calcdosestat = function(Feature, Dose_Levels, multicomp = c("none","ttest","tukey","games-howell"),
                         p.adjust.method = c("none","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),...){

# parameter check point ----------------------------------------------------------
multicomp <- match.arg(multicomp)
p.adjust.method <- match.arg(p.adjust.method)

if(sum(duplicated(Dose_Levels))>0){
  stop(paste("Duplicated Dose_Levels:",Dose_Levels[duplicated(Dose_Levels)],collapse = ", "))
}



# Initialization ----------------------------------------------------------
num_feature <- nrow(Feature)
Feature = cbind(index=seq.int(num_feature),data.table(Feature))
Normalized_Response <- data.table(index=seq.int(num_feature))
stat <- data.table(index=seq.int(num_feature))
pvalue <- data.table(index=seq.int(num_feature))
relChange <- data.table(index=seq.int(num_feature))
Dose_Levels <- as.character(Dose_Levels)
Dose_Replicates <- rep(0,length(Dose_Levels))
names(Dose_Replicates) <- Dose_Levels
SampleNames <- character()

# Statistical Analysis Part I: Basic statistics ---------------------------------------------
Response <- data.table(index=seq.int(num_feature))
mean <- matrix(nrow = num_feature, ncol = 0)
std <- matrix(nrow = num_feature, ncol = 0)
cv <- matrix(nrow = num_feature, ncol = 0)
# extract raw data by matching column names

for(DoseX in Dose_Levels) {
ResponseX <- Feature %>% dplyr::select(contains(DoseX)) # get the raw response value from X
  if(length(ResponseX)==0){
  cat("Warning message: the following dose level does not exist:", DoseX, "\n","Omit ", DoseX, "and continue...\n")
  Dose_Replicates <- Dose_Replicates[-which(Dose_Levels==DoseX)]
  Dose_Levels <- Dose_Levels[-which(Dose_Levels==DoseX)]
  } else {
  Response = cbind(Response,ResponseX)           # data clean-up
  Dose_Replicates[DoseX] <- length(ResponseX)    # number of replicates in doseX
  meanX <- apply(ResponseX, MARGIN = 1, mean) # calculate the mean in doseX
  stdX <- apply(ResponseX, MARGIN = 1, stats::sd) # calculate the std in doseX
  cvX <- stdX / meanX # calculate the coefficient of variance
  mean <- cbind(mean,meanX)
  std <- cbind(std,stdX)
  cv <- cbind(cv, cvX)
  }
}

if(length(Dose_Levels)<3){
    cat("The number of doses is less than 3. NULL is returned.\n")
    return(NULL)
}

colnames(mean) <- sapply(Dose_Levels, function(x) paste("mean",x,sep = ""))
colnames(std) <- sapply(Dose_Levels, function(x) paste("std",x,sep = ""))
colnames(cv) <- sapply(Dose_Levels, function(x) paste("cv",x, sep = ""))
SampleNames = colnames(Response[,-1L])

SampleInfo = cbind(SampleNames, dose_levels = rep(Dose_Levels, times = Dose_Replicates))
stat <- cbind(stat, mean, std, cv) # FUNCTION OUTPUT

# check if any duplicated match
if(sum(duplicated(SampleNames))>0){
stop("Following samples have duplicated matches to multiple Dose_Levels: \n",paste0(unique(SampleNames[duplicated(SampleNames)]),sep="\n"),"  Change Dose_Levels or column names for unique match.")
}

# Data normalization -------------------------------------------------
Normalized_Response <- normalization(Response[,-1L], Norm.method = "range") # "range' scaling by default
Normalized_Response <- cbind(index=Response$index, Normalized_Response)
Normalized_Response[is.na(Normalized_Response)]<-0

# Statistical Analysis Part II: Calulate pvalues and relative changes--------------------------------------------
dose_levels <- factor(rep(Dose_Levels,times=Dose_Replicates),levels=Dose_Levels)
dose_pair <- combn(unique(as.character(dose_levels)),2, paste,collapse="&")

# anova
aov_pvalue <- apply(Normalized_Response[,-1L], MARGIN =1, function(x) {a=stats::aov(as.numeric(x)~dose_levels); return(as.numeric(unlist(summary(a))["Pr(>F)1"]))})

# statictical test for p value
if(multicomp == "none"){
pairwise_pvalue <- apply(Normalized_Response[,-1L], MARGIN = 1, function(x) {
                  p <- stats::pairwise.t.test(x,dose_levels, pool.sd=FALSE, p.adjust.method="none")$p.value
                  p <- p[lower.tri(p,diag=TRUE)]; return(p)})
} else if(multicomp =="ttest"){
pairwise_pvalue <- apply(Normalized_Response[,-1L], MARGIN = 1, function(x) {
                  p <- stats::pairwise.t.test(x,dose_levels, pool.sd=FALSE, p.adjust.method=p.adjust.method)$p.value
                  p <- p[lower.tri(p,diag=TRUE)]; return(p)})
} else{
pairwise_pvalue <- apply(Normalized_Response[,-1L], MARGIN = 1, function(x) {
                  p <- userfriendlyscience::posthocTGH(x,dose_levels, method=multicomp, conf.level=0.95)$output[[ifelse(multicomp=="tukey",1,2)]]$p; return(p)}) # extract p value: first output(1) is tukey, second(2) is games-howell
}
pairwise_pvalue <- t(pairwise_pvalue)
pairwise_pvalue[is.nan(pairwise_pvalue)] <- 1 # maximum pvalue is 1.
colnames(pairwise_pvalue) <- sapply(dose_pair,function(x) paste("p_",x,sep=""))

# calculate pairwise relative changes based on normalized data
norm <- Normalized_Response[,-1L]
relchange <- matrix(nrow=nrow(norm),ncol=0)
for (i in 1:(length(Dose_Levels)-1)){
end_1 <- sum(Dose_Replicates[1:i])
start_1 <- end_1-as.numeric(Dose_Replicates[i])+1
current <- norm[,start_1:end_1]
  for(j in (i+1):length(Dose_Levels)){
  end_2 <- sum(Dose_Replicates[1:j])
  start_2 <- end_2-as.numeric(Dose_Replicates[j])+1
  higher <- norm[,start_2:end_2]
  relchange <- cbind(relchange,rowMeans(higher)-rowMeans(current))
  }
}
colnames(relchange) <- sapply(dose_pair,function(x) paste("rel_",x,sep=""))

if(multicomp=="none"){
seq_ind <- as.numeric()
j <- 1
  for(i in 1:(length(Dose_Levels)-1)){
  seq_ind <- c(seq_ind,j)
  j <- j + length(Dose_Levels)-i
  }
pairwise_pvalue<- pairwise_pvalue[,seq_ind]
relchange <- relchange[,seq_ind]
dose_pair <- dose_pair[seq_ind]
}
#rename relChange and pvalue columns

relChange <- cbind(relChange,relchange)
pvalue <- cbind(pvalue,aov_pvalue=aov_pvalue,pairwise_pvalue)
DoseStat = list(Feature = Feature, Normalized_Response = Normalized_Response, stat = stat,
                  pvalue = pvalue, relChange = relChange, Dose_Levels = Dose_Levels, Dose_Replicates = Dose_Replicates,
                  SampleInfo = SampleInfo, norm.method="range", multicomp=multicomp, p.adjust.method = p.adjust.method)

return(DoseStat)
}

#' Row-wise data normalization and transformation
#'
#' @description Perform row-wise data normalization or transformation of a given data.matrix. Return a data.matrix with normalized values.
#' @usage normalization(data,Norm.method="range",...)
#' @param data a matrix of numerical values.
#' @param Norm.method normalization methods to be applied. By default using range scaling. See details for further information.
#' @param ... Further arguments to be passed to normalization.
#' @details normalization() applies a series of commonly used metrices for data normalization. Range scaling ("range") focuses on data correlation and restricted the values in between 0 to 1, (x-min)/(max-min);
#' Auto scaling ("auto") focuses on data correlation and data is centered to 0 with a standard deviation of 1, (x-mean)/std; Pareto scaling ("pareto") focuses on data correlation and discriminates large fold-chanegs, (x-mean)/(sqrt(std));
#' Vast scaling ("vast") focuses on less varying data and discriminates largely varying data, (x-mean)/(std*cv); Level scaling ("level") focuses on fold changes against the mean value, (x-mean)/mean;
#' log10 or log2 transformation ("log10","log2") focuses on scaling the exponential relationship to a linear model, and centering data at the mean; squart root ("sqrt") is a pesudo scaling for positive values only.
#' @importFrom stats sd

normalization <- function(data,Norm.method="range",...){
  data <- as.matrix(data)
  mean <- apply(data,MARGIN = 1,mean)
  std <- apply(data,MARGIN=1,stats::sd)
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
     #cat("Notice:",paste("'",Norm.method,"'",sep=""),"method applied for normalization.\n")
  } else stop(Norm.method, " is not a normalization method See ?calcDoseStat.")
  return(norm)
}


