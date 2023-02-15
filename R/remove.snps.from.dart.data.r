#' Remove snps from a dart data object, based on a vector of locus names or locus indices
#'
#' remove.snps.from.dart.data() removes a set of nominated SNPs from a dart data object
#'
#' @param dart_data   -- dart data list               [required]
#' @param snp_indices -- indices of snps for removal  [required]
#' @param input_as_names -- if TRUE: names; if FALSE: indices
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' remove.snps.from.dart.data(dart_data, snps_to_remove, input_as_names=FALSE)
#

remove.snps.from.dart.data <- function(dart_data, snps_to_remove, input_as_names=FALSE) {

   dart_data_proc <- dart_data
   num_snps_to_remove <- length(snps_to_remove)

   if (input_as_names) {

       snp_indices <- which( colnames(dart_data_proc$gt) %in% snps_to_remove)

       if ( length(snp_indices) != num_snps_to_remove ) {
          cat(" Fatal Error: could not find all the snps that need to be remove, based on their names \n"); stop()         
       } 

   } else {
       snp_indices <- snps_to_remove
   }

   dart_data_proc$gt <- dart_data$gt[ , -snp_indices ] 
   dart_data_proc$locus_names <- dart_data$locus_names[ -snp_indices ] 
   dart_data_proc$locus_nuc <- dart_data$locus_nuc[ -snp_indices ] 
   dart_data_proc$locus_pos <- dart_data$locus_pos[ -snp_indices ] 
   dart_data_proc$locus_repro <- dart_data$locus_repro[ -snp_indices ] 

   return(dart_data_proc)

}
