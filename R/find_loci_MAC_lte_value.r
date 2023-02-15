#' Find loci that have a minor allele count (MAC) lower than (or equal to) a nominated value
#' Note: if value is set to 0, it will find fixed (non-polymorphic) loci. If none found, returns "None".
#'
#' @param dart_data       -- dart data list  [required]
#' @param value           -- minimum value of minor allele count
#' @param return_as_names -- if TRUE, returns names of loci, if false, returns indices
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' find.loci.MAC.lte.value(dart_data, value=0, input_as_names=FALSE)
#

find_loci_MAC_lte_value <- function(dart_data, value=0, return_as_names=FALSE, sample_inds) {

   # do some checking
   if(dart_data$encoding == "altcount") {
          cat(" Data found in altcount encoding, search for fixed sites commencing \n")         
   } else {
          cat(" Fatal error: MAC calculates assume altcount encoding, these calculations need to stop \n"); stop()         
   } 

   # do calculations of allele counts
   x <- dart_data$gt

   if ( !missing(sample_inds) ) {
      x <- x[ as.vector(sample_inds), ]
   }


   f0hom <- function(x) {
      length(which(x == 0))
   }

   fhet <- function(x) {
      length(which(x == 1))
   }

   f1hom <- function(x) {
      length(which(x == 2))
   }

   hom0 <- apply(x, 2, f0hom)
   het  <- apply(x, 2, fhet)
   hom1 <- apply(x, 2, f1hom)
 
   s0 <- 2*hom0 + het
   s1 <- 2*hom1 + het

   # identify loci with s0 or s1 lte value
   value0 <- which(s0 <= value)
   value1 <- which(s1 <= value)
   MAC_lte_value <- union(value0, value1)
   num_MAC_lte_value <- length(MAC_lte_value)

   # return snps with MAC lte value in appropriate format, or None
   if ( num_MAC_lte_value > 0 ) {

      if(return_as_names) {
         type <- "names"
         MAC_lte_value <- colnames(x)[MAC_lte_value]
         cat(" Found ", num_MAC_lte_value," loci with MAC lte ", value, " returning their names \n")
      } else {
         type <- "indices"
         cat(" Found ", num_MAC_lte_value," loci with MAC lte ", value, " returning their indices \n")
      }
   } else {
      cat(" Found no loci with MAC lte ", value, " returning None \n")
      type <- "none"
      MAC_lte_value <- NA
   } 

   MAC_lte_value_list <- list(MAC_lte_value=MAC_lte_value, type=type)
   
   return(MAC_lte_value_list)
}
