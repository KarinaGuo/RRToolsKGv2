#' Remove fixed SNP sites from a gentype list object  
#'
#' @param genotype list object [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' remove_fixed_snps(dart_data)
#

remove_fixed_snps <- function(dart_data) {

   dart_data_proc  <- dart_data
   genotypes       <- dart_data$gt

   fixed_snps <- find_loci_MAC_lte_value(dart_data, value=0)

   if (fixed_snps$type == "none") {
      cat("   No fixed SNPs found, no further modifications of genotype object necessary \n")
   } else {

      num_fixed <- length(fixed_snps$MAC_lte_value)
      cat("   Found ", num_fixed, "fixed SNPs, these will be removed from genotype object \n")
      dart_data_proc <- remove.snps.from.dart.data(dart_data, fixed_snps$MAC_lte_value, input_as_names=FALSE)
   }
   return(dart_data_proc)
}
