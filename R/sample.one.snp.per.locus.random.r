#' Sample one SNP per locus at random 
#'
#' sample.one.snp.per.locus.random() receives data in list format, and samples one snp per locus randomly
#' A new dart data list object is written, with the treatment field updated to "1SNPperClone"
#'
#' @param dart_data -- dart data list  [required]
#' @param seed      -- random seed    
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' sample.one.snp.per.locus.random(dart_data, seed=12345)
#
 
sample.one.snp.per.locus.random <- function(dart_data, seed=12345) {
   
   cat("   Choosing one snp per CloneID at random\n")
   cat("   Using random seed:", seed, "\n")
   set.seed(seed, kind = NULL, normal.kind = NULL)
 
   # checks #
   genotypes   <- dart_data$gt
   locus_names <- dart_data$locus_names
   num_samples <- nrow(genotypes)
   num_loci    <- ncol(genotypes)

   if( exists("genotypes") & num_samples >= 1 & num_loci >= 1) {
      cat("  QC calculations starting \n")
   } else {
      cat(" Fatal Error: could not identify a genotype table...\n"); stop()
   }

   if (length(locus_names) != num_loci) {
      cat(" Fatal Error: mismatch in size of genotype matrix, vector of locus names\n"); stop()
   } else {
      cat(" genotypes, lous names read \n")
   }

   clone_counts <- as.matrix(table(as.character(locus_names)))
   clones_with_multiple_snps <- rownames(clone_counts)[which(clone_counts > 1)]

   num_clones_with_multiple_snps <- length(clones_with_multiple_snps)
  
   if (num_clones_with_multiple_snps > 1) {
       cat(" found ", num_clones_with_multiple_snps,"clones with multiple SNPs \n")
       cat(" choosing one snp at random from each, which could take a while \n")


   for (n in 1:num_clones_with_multiple_snps) {
      
      clone             <- clones_with_multiple_snps[n]
      snps_in_clone     <- which(locus_names == clone)
      num_snps_in_clone <- length(snps_in_clone)

      if (num_snps_in_clone > 1) {
         rand_snp       <- sample(1:num_snps_in_clone)[1]
         snps_to_remove <- snps_in_clone[ -rand_snp ]
         if (n == 1) { snps_to_remove_cumulative <- snps_to_remove } 
         if (n > 1) { snps_to_remove_cumulative <- c(snps_to_remove_cumulative, snps_to_remove) } 
      } else {
         cat(" Error: this clone only contains one SNP... why? \n"); stop();
      }
   }

   cat(" A new SNP set has been selected, writing to a data object \n")

   # start building a new dart list object, to return
   dart_data_1rspl             <- dart_data
   dart_data_1rspl$treatment   <- paste(dart_data_1rspl$treatment, "1SNPperClone", sep="_")
   dart_data_1rspl$gt          <- dart_data_1rspl$gt[, -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_nuc   <- dart_data_1rspl$locus_nuc[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_names <- dart_data_1rspl$locus_names[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_pos   <- dart_data_1rspl$locus_pos[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_repro <- dart_data_1rspl$locus_repro[ -snps_to_remove_cumulative ]

   return(dart_data_1rspl)

   } else {
       cat(" there was already only one snp per clone... is this correct? \n")
       cat(" returning the data that were input with treatment appended \n")
       dart_data$treatment   <- paste(dart_data$treatment, "1SNPperClone", sep="_")
       return(dart_data)
   }


}
