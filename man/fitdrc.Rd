% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitdrc.R
\name{fitdrc}
\alias{fitdrc}
\title{Monotonic trend fitting and ED50 estimation for metabolomic dose-reponse analysis (TOXcms step 3)}
\usage{
fitdrc(DoseResponse_report, Dose_values, ED=0.5, mz_tag = "mzmed", rt_tag = "rtmed",
export = TRUE, plot = TRUE,...)
}
\arguments{
\item{DoseResponse_report}{list The output list object from trendfilter().}

\item{Dose_values}{numeric A vector of dose values.}

\item{ED}{numeric A value between (0,1), indicating the effective dose.}

\item{mz_tag}{character Name of the m/z column in the feature table.}

\item{rt_tag}{character Name of the retention time column in the feature table.}

\item{export}{TRUE/FALSE. When set to TRUE, a csv file including all results will be exported.}

\item{plot}{TRUE/FALSE. When set to TRUE, a PDF file including all plots will be exported.}

\item{...}{Further arguments to be passed to fitdrc}
}
\value{
The function returns doseresponse_report, a list object including same elements as the input and an addition of ED50 value.
}
\description{
The fitdrc() fits the monotonic metabolic trends to a 4-parameter logistic model and visualize the fitting curve. Based on the fitted model,
    the function calculates the ED50 value for metabolic trend. The function returns the same list object as the input with additional table of ED50 values.
    The function also exports a PDF file including all plots. The output object is called doseResponse_report.
}
\author{
Cong-Hui Yao <conghui.yao@wustl.edu>
    Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
}
