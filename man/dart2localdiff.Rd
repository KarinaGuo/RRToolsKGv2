% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dart2localdiff.r
\name{dart2localdiff}
\alias{dart2localdiff}
\title{Prepares input files that are ready
for analysis using software localdiff}
\usage{
dart2localdiff(dms, basedir, species, dataset, pop)
}
\arguments{
\item{dms}{-- a dart data object [Required]}

\item{basedir}{-- name of the base directory for R&R}

\item{species}{-- species name}

\item{dataset}{-- dataset name}

\item{pop}{-- a vector containing population assignments for each sample}
}
\value{
a list containing names of infiles for localdiff analysis
}
\description{
Input is a dart data object with altcount genotype encoding
(0=hom ref, 1=het, 2=hom alt).
}
\examples{
\dontrun{
ld_files <- dart2localdiff(gms, etc)
}
}
\author{
Jason Bragg (jasongbragg@gmail.com)
}
