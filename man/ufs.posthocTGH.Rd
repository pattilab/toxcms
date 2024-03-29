% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{ufs.posthocTGH}
\alias{ufs.posthocTGH}
\title{Tukey and Games-Howell post-hoc analysis adopted from userfriendlyscience}
\usage{
ufs.posthocTGH(y, x, method = c("games-howell", "tukey"), 
conf.level = 0.95, digits = 2, p.adjust = "none", formatPvalue = TRUE)
}
\arguments{
\item{y}{numeric a numeric vector to be tested. In toxcms, this is a series of response values.}

\item{x}{factor a factor indicating sample groups. In toxcms, this is the dose levels.}

\item{method}{character post-hoc methods to be applied ("games-howell" or “tukey").}

\item{conf.level}{confidence level of the confidence intervals.}

\item{digits}{numeric(1) the number of digits to show in the output}

\item{p.adjust}{character Project name of the TOXcms analysis。 See also \link[stats]{p.adjust}}

\item{formatPvalue}{TRUE/FALSE whether to print out p values according to APA standards.}
}
\value{
calcdosestat returns a list object consisting all the statistical results described above.
}
\description{
The ufs.posthocTGH function is adopted from userfriendlyscience package. This function is used along with 
calcdosestat function to perform post-hoc analysis using Tukey or Games-Howell approaches.
}
\author{
Cong-Hui Yao <conghui.yao@wustl.edu>
Lingjue Wang (Mike) <wang.lingjue@wustl.edu>
}
