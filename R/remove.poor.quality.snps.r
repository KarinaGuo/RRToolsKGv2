#' Removes SNP sites from a DArT data object subject to quality criteria  
#'
#' remove.poor.quality.snps() receives dart data in list format, 
#' and removes SNP sites that fall below presence and RepAvg thresholds
#'
#' @param dart_data -- dart data list  [required]
#' @param min_repro   -- minimum value of RepAvg
#' @param max_missing -- maximum proportion of missing data [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' remove.poor.quality.snps(dart_data, min_repro=0.975, max_missing=0.25)
#

remove.poor.quality.snps <- function(dart_data, min_repro=0.975, max_missing=0.25) {

   bad_snps_repro   <- find.loci.repro.lte.value(dart_data, value=min_repro, return_as_names=FALSE)
   bad_snps_missing <- find.loci.missing.gte.value(dart_data, value=max_missing, return_as_names=FALSE)

   if (bad_snps_repro$type == "none" & bad_snps_missing$type == "none") {

      cat(" No bad SNPs found, no further modifications of DArT object necessary \n")
      return(dart_data)

   } else if (bad_snps_repro$type == "indices" | bad_snps_missing$type == "indices") {

      if (bad_snps_repro$type == "indices" & bad_snps_missing$type == "indices") {
         cat(" Returning missing and low reproducibility snps as bad \n")
         bad_snp_indices <- union(bad_snps_repro$repro_lte_value, bad_snps_missing$missing_gte_value)
      } 

      if (bad_snps_repro$type == "none" & bad_snps_missing$type == "indices") {
         cat(" Returning missing snps as bad \n")
         bad_snp_indices <- bad_snps_missing$missing_gte_value
      } 

      if (bad_snps_repro$type == "indices" & bad_snps_missing$type == "none") {
         cat(" Returning low reproducibility snps as bad \n")
         bad_snp_indices <- bad_snps_repro$repro_lte_value
      } 

      dart_data_proc <- remove.snps.from.dart.data(dart_data, bad_snp_indices, input_as_names=FALSE)
      dart_data_proc$treatment   <- paste(dart_data$treatment, "SNPFilt", sep="_")
      return(dart_data_proc)

   } else {
      cat(" Fatal error: did not return useable quality info... \n"); stop() 
   }

}
