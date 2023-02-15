#' Sample one SNP per locus at random 
#'
#' downsample.snps() receives data in list format, and samples a nominated number
#' of snp sites at random, with the treatment field updated to "downSamp"
#'
#' @param dart_data -- dart data list  [required]
#' @param number    -- number of random snps to select
#' @param seed      -- random seed    
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' downsample.snps(dart_data, number=1000, seed=12345)
#

downsample_snps <- function(dart_data, number=1000, seed=12345) {

   set.seed(seed, kind = NULL, normal.kind = NULL)
   dart_data_proc   <- dart_data
   num_snps_in_data <- ncol(dart_data$gt)

   # checks
   if (number >= num_snps_in_data) {
      cat(" Error: the number of target snps is fewer or equal to number available \n");  
      cat(" Returning the data, without modification \n");  
      dart_downsampled <- dart_data_proc
   } else {
      cat("   Choosing ", number, "snp sites at random\n")
      cat("   Using random seed:", seed, "\n")
      snps_to_remove <- sample(1:num_snps_in_data)[-c(1:number)]      
      dart_downsampled <- remove.snps.from.dart.data(dart_data_proc, snps_to_remove, input_as_names=FALSE)
   }
   dart_downsampled$treatment   <- paste(dart_data_proc$treatment, "downSamp", sep="_")
   return(dart_downsampled)

}
