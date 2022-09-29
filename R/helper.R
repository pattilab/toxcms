#' Tukey and Games-Howell post-hoc analysis adopted from userfriendlyscience
#'
#' @description The ufs.posthocTGH function is adopted from userfriendlyscience package. This function is used along with 
#' calcdosestat function to perform post-hoc analysis using Tukey or Games-Howell approaches.
#' @usage ufs.posthocTGH(y, x, method = c("games-howell", "tukey"), 
#' conf.level = 0.95, digits = 2, p.adjust = "none", formatPvalue = TRUE)
#' @param y numeric a numeric vector to be tested. In toxcms, this is a series of response values.
#' @param x factor a factor indicating sample groups. In toxcms, this is the dose levels. 
#' @param method character post-hoc methods to be applied ("games-howell" or “tukey").
#' @param conf.level confidence level of the confidence intervals.
#' @param digits numeric(1) the number of digits to show in the output
#' @param p.adjust character Project name of the TOXcms analysis。 See also \link[stats]{p.adjust}
#' @param formatPvalue TRUE/FALSE whether to print out p values according to APA standards.
#' @return calcdosestat returns a list object consisting all the statistical results described above.
#' @author Cong-Hui Yao <conghui.yao@wustl.edu>
#' Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
#' @import utils stats dplyr
#' @export

ufs.posthocTGH <- function (y, x, method = c("games-howell", "tukey"), 
          conf.level = 0.95, digits = 2, p.adjust = "none", formatPvalue = TRUE) 
{
  method <- tolower(method)
  tryCatch(method <- match.arg(method), error = function(err) {
    stop("Argument for 'method' not valid!")
  })
  res <- list(input = as.list(environment()))
  res$intermediate <- list(x = factor(x[complete.cases(x, y)]), 
                           y = y[complete.cases(x, y)])
  res$intermediate$n <- tapply(y, x, length)
  res$intermediate$groups <- length(res$intermediate$n)
  res$intermediate$df <- sum(res$intermediate$n) - res$intermediate$groups
  res$intermediate$means <- tapply(y, x, mean)
  res$intermediate$variances <- tapply(y, x, var)
  res$intermediate$names <- levels(res$intermediate$x)
  res$intermediate$pairNames <- combn(res$intermediate$groups, 
                                      2, function(ij) {
                                        paste0(rev(res$intermediate$names[ij]), collapse = "-")
                                      })
  res$intermediate$descriptives <- cbind(res$intermediate$n, 
                                         res$intermediate$means, res$intermediate$variances)
  rownames(res$intermediate$descriptives) <- levels(res$intermediate$x)
  colnames(res$intermediate$descriptives) <- c("n", "means", 
                                               "variances")
  res$intermediate$errorVariance <- sum((res$intermediate$n - 
                                           1) * res$intermediate$variances)/res$intermediate$df
  res$intermediate$se <- combn(res$intermediate$groups, 2, 
                               function(ij) {
                                 sqrt(res$intermediate$errorVariance * sum(1/res$intermediate$n[ij]))
                               })
  res$intermediate$dmeans <- combn(res$intermediate$groups, 
                                   2, function(ij) {
                                     diff(res$intermediate$means[ij])
                                   })
  res$intermediate$t <- abs(res$intermediate$dmeans)/res$intermediate$se
  res$intermediate$p.tukey <- ptukey(res$intermediate$t * sqrt(2), 
                                     res$intermediate$groups, res$intermediate$df, lower.tail = FALSE)
  res$intermediate$alpha <- (1 - conf.level)
  res$intermediate$qcrit <- qtukey(res$intermediate$alpha, 
                                   res$intermediate$groups, res$intermediate$df, lower.tail = FALSE)/sqrt(2)
  res$intermediate$tukey.low <- res$intermediate$dmeans - (res$intermediate$qcrit * 
                                                             res$intermediate$se)
  res$intermediate$tukey.high <- res$intermediate$dmeans + 
    (res$intermediate$qcrit * res$intermediate$se)
  res$output <- list()
  res$output$tukey <- data.frame(res$intermediate$dmeans, res$intermediate$tukey.low, 
                                 res$intermediate$tukey.high, res$intermediate$t, res$intermediate$df, 
                                 res$intermediate$p.tukey)
  columnNames <- c("diff", "ci.lo", "ci.hi", 
                   "t", "df", "p")
  if (p.adjust != "none") {
    res$output$tukey$p.tukey.adjusted <- p.adjust(res$intermediate$p.tukey, 
                                                  method = p.adjust)
    columnNames <- c(columnNames, "p.adjusted")
  }
  rownames(res$output$tukey) <- res$intermediate$pairNames
  colnames(res$output$tukey) <- columnNames
  res$intermediate$df.corrected <- combn(res$intermediate$groups, 
                                         2, function(ij) {
                                           sum(res$intermediate$variances[ij]/res$intermediate$n[ij])^2/sum((res$intermediate$variances[ij]/res$intermediate$n[ij])^2/(res$intermediate$n[ij] - 
                                                                                                                                                                         1))
                                         })
  res$intermediate$se.corrected <- combn(res$intermediate$groups, 
                                         2, function(ij) {
                                           sqrt(sum(res$intermediate$variances[ij]/res$intermediate$n[ij]))
                                         })
  res$intermediate$t.corrected <- abs(res$intermediate$dmeans)/res$intermediate$se.corrected
  res$intermediate$qcrit.corrected <- qtukey(res$intermediate$alpha, 
                                             res$intermediate$groups, res$intermediate$df.corrected, 
                                             lower.tail = FALSE)/sqrt(2)
  res$intermediate$gh.low <- res$intermediate$dmeans - res$intermediate$qcrit.corrected * 
    res$intermediate$se.corrected
  res$intermediate$gh.high <- res$intermediate$dmeans + res$intermediate$qcrit.corrected * 
    res$intermediate$se.corrected
  res$intermediate$p.gameshowell <- ptukey(res$intermediate$t.corrected * 
                                             sqrt(2), res$intermediate$groups, res$intermediate$df.corrected, 
                                           lower.tail = FALSE)
  res$output$games.howell <- data.frame(res$intermediate$dmeans, 
                                        res$intermediate$gh.low, res$intermediate$gh.high, res$intermediate$t.corrected, 
                                        res$intermediate$df.corrected, res$intermediate$p.gameshowell)
  columnNames <- c("diff", "ci.lo", "ci.hi", 
                   "t", "df", "p")
  if (p.adjust != "none") {
    res$output$games.howell$p.gameshowell.adjusted <- p.adjust(res$intermediate$p.gameshowell, 
                                                               method = p.adjust)
    columnNames <- c(columnNames, "p.adjusted")
  }
  rownames(res$output$games.howell) <- res$intermediate$pairNames
  colnames(res$output$games.howell) <- columnNames
  class(res) <- "posthocTGH"
  return(res)
}

#' Row-wise data normalization and transformation
#'
#' @description Perform row-wise data normalization or transformation of a given data.matrix. Return a data.matrix with normalized values.
#' @usage normalization(data,Norm.method=c("range","auto","pareto","vast","level",
#' "log10","log2","sqrt"),...)
#' @param data a matrix of numerical values.
#' @param Norm.method normalization methods to be applied. By default using range scaling. See details for further information.
#' @param ... Further arguments to be passed to normalization.
#' @details normalization() applies a series of commonly used metrices for data normalization. Range scaling ("range") focuses on data correlation and restricted the values in between 0 to 1, (x-min)/(max-min);
#' Auto scaling ("auto") focuses on data correlation and data is centered to 0 with a standard deviation of 1, (x-mean)/std; Pareto scaling ("pareto") focuses on data correlation and discriminates large fold-chanegs, (x-mean)/(sqrt(std));
#' Vast scaling ("vast") focuses on less varying data and discriminates largely varying data, (x-mean)/(std*cv); Level scaling ("level") focuses on fold changes against the mean value, (x-mean)/mean;
#' log10 or log2 transformation ("log10","log2") focuses on scaling the exponential relationship to a linear model, and centering data at the mean; squart root ("sqrt") is a pesudo scaling for positive values only.
#' @importFrom stats sd
#' @import magrittr

normalization <- function(data,Norm.method=c("range","auto","pareto","vast","level","log10","log2","sqrt"),...){

  Norm.method <- match.arg(Norm.method)
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

#' Calculate observed ED values using dr4pl-fitted parameters.
#' 
#' @usage ObservedED(ED, theta)
#' @param ED numeric the level of effective dose from 0 to 1.
#' @param theta the parameter set of logistic regression curve.

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
