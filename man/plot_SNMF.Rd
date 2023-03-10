% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_SNMF.r
\name{plot_SNMF}
\alias{plot_SNMF}
\title{Make a plot of admixture coefficients from an sNMF analysis}
\usage{
plot_SNMF(basedir, species, dataset, treatment, Kvals = 2:6, Rvals = 1:3)
}
\arguments{
\item{basedir}{-- Base dirctory   [required]}

\item{species}{-- species name    [required]}

\item{dataset}{-- dataset identifier  [required]}

\item{treatment}{-- treatment       [required]}

\item{K}{-- sNMF K value    [required]}

\item{r}{-- sNMF replicate  [required]}
}
\description{
plot_SNMF() information about an sNMF run, opens the project
A pdf file of admixture coefficients is made for a nominated K, R
}
\author{
Jason Bragg (jasongbragg@gmail.com)
}
