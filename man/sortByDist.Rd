% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusttrend.R
\name{sortByDist}
\alias{sortByDist}
\title{Sorting dose-response trends by absolute distances in hierarchical cllustering.}
\usage{
sortByDist(DoseResponse_report, Index, Sort.method=c("near","far","both"),
Sort.thres = 20, Dist.method = "euclidean", Hclust.method = "average",...)
}
\arguments{
\item{DoseResponse_report}{list This is the output of trendfilter(), the end of toxcms step 2.}

\item{Index}{numeric feature indexes to be queried.}

\item{Sort.method}{character methods to sort distances, can be one of "far","near" or "both".}

\item{Sort.thres}{numeric sorting threshold for each group based on Sort.method.}

\item{Dist.method}{Matrice of similarity calculation. By default "euclidean". See also \link[stats]{dist}.}

\item{Hclust.method}{clustering approach. By default is "average". See also \link[stats]{hclust}.}

\item{...}{Further arguments to be passed.}
}
\value{
list Each element includes a vector of reference indexes corresponding to the reference_index.
}
\description{
This function sorts similar trends of refence trend(s) based on similarity matrices.
}
\author{
Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
}
