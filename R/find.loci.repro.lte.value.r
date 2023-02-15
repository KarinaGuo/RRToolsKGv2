#' Find loci that have a reproducibility value lower than (or equal to) a nominated value
#'
#' find.loci.repro.lte.value() searches a dart data object for SNPs with reproducibility less than or equal to a nominated value
#'
#' @param dart_data       -- dart data list  [required]
#' @param value           -- minimum value of reproducibility
#' @param return_as_names -- if TRUE, returns names of loci, if false, returns indices
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' find.loci.repro.lte.value(dart_data, value=0, input_as_names=FALSE)
#

find.loci.repro.lte.value <- function(dart_data, value, return_as_names=FALSE) {

   repro <- dart_data$locus_repro

   # do some checking
   if( length(repro) > 1 & length(repro) == ncol(dart_data$gt) ) {
          cat(" Data reproducibility stats found, finding values less than ", value, "\n")         
   } else {
          cat(" Fatal error: did not find reproducibility values \n"); stop()         
   } 

   repro_lte_value <- which(repro <= value)
   num_repro_lte_value <- length(repro_lte_value)

   x <- dart_data$gt

   # return snps with MAC lte value in appropriate format, or None
   if ( num_repro_lte_value > 0 ) {

      if(return_as_names) {
         type <- "names"
         repro_lte_value <- colnames(x)[repro_lte_value]
         cat(" Found ", num_repro_lte_value," loci with reproducibility <= ", value, " returning their names \n")
      } else {
         type <- "indices"
         cat(" Found ", num_repro_lte_value," loci with reproducibility <= ", value, " returning their indices \n")
      }
   } else {
      cat(" Found no loci with reproducibility <= ", value, " returning None \n")
      type <- "none"
      repro_lte_value <- NA
   } 

   repro_lte_value_list <- list(repro_lte_value=repro_lte_value, type=type)
   
   return(repro_lte_value_list)
}
