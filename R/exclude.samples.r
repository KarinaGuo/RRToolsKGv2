#' Remove a list of samples from a DArT data object  
#' Optionally removes markers that become monomorphic
#'
#' @param dart_data              -- dart data list  [required]
#' @param excluded_sample_file   -- name of a file containing list of samples to be excluded 
#' @param removed_fixed_loci     -- removed snps that are fixed after removal of nominated samples
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' exclude.samples(dart_data, excluded_sample_file, remove_fixed_loci=TRUE)
#

exclude.samples <- function(dart_data, by="file", excluded_sample_file="None", excluded_sample_names, remove_fixed_loci=TRUE) {

   dart_data_proc <- dart_data
   genotypes      <- dart_data$gt

   if ( by == "file" | by == "names" ) {

      if (by == "file") {
         samples_to_exclude <- as.matrix(read.table(excluded_sample_file, header=FALSE))
      } 
 
      if (by == "names") {
         samples_to_exclude <- excluded_sample_names
      } 
      

      if ( exists("samples_to_exclude") ) {

         num_excluded <- length(samples_to_exclude)
         cat("  Proceeding with exclusion of ", num_excluded, " samples listed in ", excluded_sample_file, "\n")
         ind_remove <- which( rownames(genotypes) %in% samples_to_exclude)

         # check that all nominated for removal were found
         # if so, remove and create an updated version of dart_data
         # if not, issue warning ... 
         if ( length(ind_remove) == num_excluded ) {

              dart_data_proc$gt           <- genotypes[ -ind_remove, ]
              dart_data_proc$sample_names <- rownames(genotypes)[ -ind_remove ]
              dart_data_proc$treatment    <- paste(dart_data_proc$treatment, "sampFilt", sep="_")

         } else {
            cat("  Warning: Could not find all the samples that were nominated for exclusion in the data matrix -- check names \n")    
         }

         if (remove_fixed_loci) {
            cat("  Proceeding with removal of loci that become fixed \n")
            dart_data_proc <- remove.fixed.snps(dart_data_proc)
             
         } else {
            cat("  Warning: samples were removed, and some loci might be fixed in the remaining samples. Consider removing non-polymorphic sites \n")
         }

      } else {
         cat("  Fatal error: could not read sample names from ", excluded_sample_file, "\n"); stop();
      }


   } else {
      cat("  No samples need to be excluded based on missingness of genotypes \n ")
   } # if excluded_sample_file equals "None"

   return(dart_data_proc)
}
