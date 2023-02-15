#' Returns loci that are near each other on a eucalypt chromosome
#' 
#' Input is a dart, meta, genome object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a list containing info about snp pairs
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

get.intra.locus.pairs <- function(dart_data, MAF_thresh=0.2, min_dist_thresh_bp=10) {
   
   AF_summary <- calculate.AF.summary(dart_data)

   locus_names <- dart_data$locus_names
   locus_pos   <- dart_data$locus_pos
   genotypes   <- dart_data$gt

   clone_counts         <- as.matrix(table(as.character(locus_names)))
   clones_with_two_snps <- rownames(clone_counts)[which(clone_counts == 2)]

   num_clones_with_two_snps <- length(clones_with_two_snps)

   count <- 1
   for (c in 1:num_clones_with_two_snps) {

      clone <- clones_with_two_snps[c]
      ind_clone <- which(locus_names == clone)
      snp_dist  <- abs(locus_pos[ ind_clone[1] ] - locus_pos[ ind_clone[2] ])

      if( AF_summary$P[ ind_clone[1] ] > MAF_thresh & AF_summary$P[ ind_clone[1] ] < (1-MAF_thresh) & AF_summary$P[ ind_clone[2] ] > MAF_thresh & AF_summary$P[ ind_clone[2] ] < (1-MAF_thresh) ) {
         pair <- c( colnames(genotypes)[ind_clone[1]], colnames(genotypes)[ind_clone[2]], snp_dist)
         if (count == 1) {pairs <- pair}
         if (count > 1) {pairs <- rbind(pairs, pair)}
         count <- count + 1
      }

   }

   pairs <- pairs[which( as.numeric(pairs[,3]) > min_dist_thresh_bp),]
   return(pairs)
}



