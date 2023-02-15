#' Returns loci that are near each other on a eucalypt chromosome
#' 
#' Input is a dart, meta, genome object
#' 
#' @param dart_data   -- a dart meta data object [Required]
#' MAF_thresh         -- minimum minor allele frequency
#' min_dist_thresh_bp -- minimum distance apart of loci
#' max_dist_thresh_bp -- maximum distance apart of loci
#' @return a list containing info about snp pairs
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

get.genome.locus.pairs <- function(dart_data, MAF_thresh, min_dist_thresh_bp=150, max_dist_thresh_bp=1E7) {
   
   AF_summary <- calculate.AF.summary(dart_data)
   scaffold_number <- 11
   for (s in 1:scaffold_number) {

   scaffold      <- paste("scaffold_", s, sep="")

   ind_scaff     <- which( dart_data$eucagen_locations$Chrom_Eucalyptus_v10_phytozome == scaffold )
   ind_MAF       <- which( AF_summary$P > MAF_thresh & AF_summary$P < (1-MAF_thresh) )

   ind_scaff_MAF <- intersect(ind_scaff, ind_MAF)

   dist_bp       <- dist( dart_data$eucagen_locations$ChromPos_Eucalyptus_v10_phytozome[ ind_scaff_MAF ]  )

   m_dist_bp     <- as.matrix(dist_bp)
   colnames(m_dist_bp) <- colnames(dart_data$gt)[ ind_scaff_MAF ]
   rownames(m_dist_bp) <- colnames(dart_data$gt)[ ind_scaff_MAF ]

   c <- 1
   for (i in 1:length(ind_scaff_MAF)) {
      for (j in 1:length(ind_scaff_MAF)) {

         if (j > i) { 
             d_bp <- m_dist_bp[i,j]
             if (d_bp <= max_dist_thresh_bp) {
                 pair <- c( colnames(m_dist_bp)[i], colnames(m_dist_bp)[j], d_bp)
                 if (c == 1) {pairs <- pair}
                 if (c > 1) {pairs <- rbind(pairs, pair)}
                 c <- c + 1
             }
         }

      }
   }

   }
   pairs <- pairs[which( as.numeric(pairs[,3]) > min_dist_thresh_bp),]
   return(pairs)
}
