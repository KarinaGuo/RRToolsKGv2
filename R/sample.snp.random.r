#' Sample one SNP per locus at random 
#'
#' sample.snp.random() receives data in list format, and samples a nominated number
#' of snp sites at random, with the treatment field updated to "downSamp"
#'
#' @param dart_data -- dart data list  [required]
#' @param number    -- number of random snps to select
#' @param seed      -- random seed    
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' sample.snp.random(dart_data, seed=12345)
#
 
sample.snp.random <- function(dart_data, number=5000, seed=12345) {

   seed <- 12345
   set.seed(seed, kind = NULL, normal.kind = NULL)
 
   # checks #
   genotypes   <- dart_data$gt
   locus_names <- dart_data$locus_names
   num_samples <- nrow(genotypes)
   num_loci    <- ncol(genotypes)

   if( exists("genotypes") & num_samples >= 1 & num_loci >= 1) {
      cat("  Found genotype matrix \n")
   } else {
      cat(" Fatal Error: could not identify a genotype table...\n"); stop()
   }

   if (length(locus_names) != num_loci) {
      cat(" Fatal Error: mismatch in size of genotype matrix, vector of locus names\n"); stop()
   } else {
      cat(" genotypes, lous names read \n")
   }

   if (num_loci <= number) {
      cat(" Fatal Error: number of requested SNPs is smaller than the number available \n"); stop()
   } else {
      cat("   Choosing ", number, "snp sites at random\n")
      cat("   Using random seed:", seed, "\n")
   }


   snps_to_keep   <- sample(1:num_loci)[1:number]

   cat(" A new SNP set has been selected, writing to a data object \n")

   # start building a new dart list object, to return
   dart_data_proc             <- dart_data
   dart_data_proc$treatment   <- paste(dart_data_proc$treatment, "downSamp", sep="_")
   dart_data_proc$gt          <- dart_data_proc$gt[, snps_to_keep ]
   dart_data_proc$locus_nuc   <- dart_data_proc$locus_nuc[ snps_to_keep ]
   dart_data_proc$locus_names <- dart_data_proc$locus_names[ snps_to_keep ]
   dart_data_proc$locus_pos   <- dart_data_proc$locus_pos[ snps_to_keep ]
   dart_data_proc$locus_repro <- dart_data_proc$locus_repro[ snps_to_keep ]

   return(dart_data_proc)

}
