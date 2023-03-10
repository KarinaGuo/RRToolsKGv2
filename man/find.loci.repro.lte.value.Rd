% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find.loci.repro.lte.value.r
\name{find.loci.repro.lte.value}
\alias{find.loci.repro.lte.value}
\title{Find loci that have a reproducibility value lower than (or equal to) a nominated value}
\usage{
find.loci.repro.lte.value(dart_data, value, return_as_names = FALSE)
}
\arguments{
\item{dart_data}{-- dart data list  [required]}

\item{value}{-- minimum value of reproducibility}

\item{return_as_names}{-- if TRUE, returns names of loci, if false, returns indices}
}
\description{
find.loci.repro.lte.value() searches a dart data object for SNPs with reproducibility less than or equal to a nominated value
}
\examples{
find.loci.repro.lte.value(dart_data, value=0, input_as_names=FALSE)
}
\author{
Jason Bragg (jasongbragg@gmail.com)
}
