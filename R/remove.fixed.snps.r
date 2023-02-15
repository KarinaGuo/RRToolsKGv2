#' Remove fixed SNP sites from a DArT data object  
#'
#' remove.fixed.snps() receives dart data in list format, and removes SNP sites that are fixed
#'
#' @param dart_data -- dart data list  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' remove.fixed.snps(dart_data)
#

remove.fixed.snps <- function(dart_data) {

   dart_data_proc  <- dart_data
   genotypes       <- dart_data$gt

   fixed_snps <- find.loci.MAC.lte.value(dart_data, value=0)

   if (fixed_snps$type == "none") {
      cat(" No fixed SNPs found, no further modifications of DArT object necessary \n")
   } else {

      num_fixed <- length(fixed_snps$MAC_lte_value)
      cat(" Found ", num_fixed, "fixed SNPs, these will be removed from DArT object \n")
      dart_data_proc <- remove.snps.from.dart.data(dart_data, fixed_snps$MAC_lte_value, input_as_names=FALSE)
   }
   return(dart_data_proc)
}
