% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove.snps.from.dart.data.r
\name{remove.snps.from.dart.data}
\alias{remove.snps.from.dart.data}
\title{Remove snps from a dart data object, based on a vector of locus names or locus indices}
\usage{
remove.snps.from.dart.data(dart_data, snps_to_remove, input_as_names = FALSE)
}
\arguments{
\item{dart_data}{-- dart data list               [required]}

\item{input_as_names}{-- if TRUE: names; if FALSE: indices}

\item{snp_indices}{-- indices of snps for removal  [required]}
}
\description{
remove.snps.from.dart.data() removes a set of nominated SNPs from a dart data object
}
\examples{
remove.snps.from.dart.data(dart_data, snps_to_remove, input_as_names=FALSE)
}
\author{
Jason Bragg (jasongbragg@gmail.com)
}
