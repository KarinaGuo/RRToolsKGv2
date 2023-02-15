#' Find loci that have a reproducibility value lower than (or equal to) a nominated value
#'
#' find.loci.missing.gte.value() searches a dart data object for SNPs with reproducibility less than or equal to a nominated value
#'
#' @param dart_data       -- dart data list  [required]
#' @param value           -- maximum peoportion of missing data for a snp
#' @param return_as_names -- if TRUE, returns names of loci, if false, returns indices
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' find.loci.missing.gte.value(dart_data, value=0, input_as_names=FALSE)
#

find.loci.missing.gte.value <- function(dart_data, value, return_as_names=FALSE) {

    x <- dart_data$gt
   mm <- dart_data$gt

   mm[ which(is.na(x)) ]  <- 1
   mm[ which(!is.na(x)) ] <- 0

   num_samples <- nrow(x)
   prop_missing <- colSums(mm) / num_samples
   
   # do some checking
   if( length(prop_missing) > 1 & length(prop_missing) == ncol(dart_data$gt) ) {
          cat(" Missing data stats calculated, finding values less than ", value, "\n")         
   } else {
          cat(" Fatal error: could not calculate proportion missing values \n"); stop()         
   } 

   missing_gte_value <- which(prop_missing >= value) 
   num_missing_gte_value <- length(missing_gte_value)

   # return snps with MAC lte value in appropriate format, or None
   if ( num_missing_gte_value > 0 ) {

      if(return_as_names) {
         type <- "names"
         missing_gte_value <- colnames(x)[missing_gte_value]
         cat(" Found ", num_missing_gte_value," loci with missingness >= ", value, " returning their names \n")
      } else {
         type <- "indices"
         cat(" Found ", num_missing_gte_value," loci with missingness >= ", value, " returning their indices \n")
      }
   } else {
      cat(" Found no loci with missingness >= ", value, " returning None \n")
      type <- "none"
      missing_gte_value <- NA
   } 

   missing_gte_value_list <- list(missing_gte_value=missing_gte_value, type=type)
   
   return(missing_gte_value_list)
}
