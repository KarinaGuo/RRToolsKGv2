#' Sample one SNP per cluster at random 
#'
#' @param dart_data   -- dart data list  [required]
#' @param cluster_fil -- name of cluster file [required]
#' @param seed        -- random seed    
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)

 
sample_one_snp_per_cluster_random <- function(dart_data, cluster_fil, seed=12345) {

   
   cat("   Choosing one snp per cluster at random\n")
   cat("   Using random seed:", seed, "\n")
   set.seed(seed, kind = NULL, normal.kind = NULL)
 
   # checks #
   genotypes   <- dart_data$gt
   locus_names <- colnames(dart_data$gt)
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

   clusters <- read.table(cluster_fil, header=FALSE, sep=",")
   colnames(clusters) <- c("snp", "clust")
   i_snp_clust <- match(  locus_names,  clusters$snp )

   # clusters ordered by locus names in gt data:
   snp_clust   <- clusters$clust[i_snp_clust]

   cluster_counts <- as.matrix(table(as.character(snp_clust)))
   clusters_with_multiple_snps <- rownames(cluster_counts)[which(cluster_counts > 1)]

   num_clusters_with_multiple_snps <- length(clusters_with_multiple_snps)
  
   if (num_clusters_with_multiple_snps > 1) {
       cat(" found ", num_clusters_with_multiple_snps,"clusters with multiple SNPs \n")
       cat(" choosing one snp at random from each, which could take a while \n")
   } else {
       cat(" there was already only one snp per cluster... is this correct? \n")
   }


   for (n in 1:num_clusters_with_multiple_snps) {
      
      clust             <- clusters_with_multiple_snps[n]
      snps_in_cluster   <- which(snp_clust == clust)
      num_snps_in_cluster <- length(snps_in_cluster)

      if (num_snps_in_cluster > 1) {
         rand_snp       <- sample(1:num_snps_in_cluster)[1]
         snps_to_remove <- snps_in_cluster[ -rand_snp ]
         if (n == 1) { snps_to_remove_cumulative <- snps_to_remove } 
         if (n > 1) { snps_to_remove_cumulative <- c(snps_to_remove_cumulative, snps_to_remove) } 
      } else {
         cat(" Error: this cluster only contains one SNP... why? \n"); stop();
      }
   }

   cat(" A new SNP set has been selected, writing to a data object \n")

   # start building a new dart list object, to return
   dart_data_1rspl             <- dart_data
   dart_data_1rspl$treatment   <- paste(dart_data_1rspl$treatment, "1SNPperCluster", sep="_")
   dart_data_1rspl$gt          <- dart_data_1rspl$gt[, -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_nuc   <- dart_data_1rspl$locus_nuc[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_names <- dart_data_1rspl$locus_names[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_pos   <- dart_data_1rspl$locus_pos[ -snps_to_remove_cumulative ]
   dart_data_1rspl$locus_repro <- dart_data_1rspl$locus_repro[ -snps_to_remove_cumulative ]

   return(dart_data_1rspl)

}
