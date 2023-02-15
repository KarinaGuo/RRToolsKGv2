#' Returns estimates of Linkage
#' 
#' Input is a dart, meta, genome object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a list containing info about snp pairs
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

calculate.genome.LD.decay <- function(dart_data, MAF_thresh=0.1, min_dist_thresh_bp=150, max_dist_thresh_bp=1E7) {

   # load some required libraries
   require(genetics)
   require(LDcorSV)

   #AF_summary <- calculate.AF.summary(dart_data)
   #MAF_thresh <- 0.2
   #dist_thresh_bp  <- 1E7

   pairs <- get.genome.locus.pairs(dart_data, MAF_thresh, min_dist_thresh_bp, max_dist_thresh_bp) 
   npairs <- nrow(pairs)

   altcount2genotype <- function(vector) {

      v <- vector
      v[ vector == 0 ] <- "R/R"
      v[ vector == 1 ] <- "R/A"
      v[ vector == 2 ] <- "A/A"
     
      g <- as.genotype(v)
      return(g)
   }


   for (p in 1:npairs) {

      i_L1 <- which( colnames(dart_data$gt) == pairs[p,1]) 
      i_L2 <- which( colnames(dart_data$gt) == pairs[p,2]) 

      v1 <- dart_data$gt[,i_L1]
      v2 <- dart_data$gt[,i_L2]
 
      g1 <- altcount2genotype( v1 )
      g2 <- altcount2genotype( v2 )

      p_LD <- LD(g1,g2)

      if (p == 1) { D_pairs <- p_LD$D }
      if (p  > 1) { D_pairs <- c(D_pairs, p_LD$D) }

      LDcorSV_object <- dart_data$gt[, c(i_L1, i_L2)]
      R2 <- Measure.R2(LDcorSV_object, na.presence=TRUE)

      if (p == 1) { R2_pairs <- R2 }
      if (p  > 1) { R2_pairs <- c(R2_pairs, R2) }

   }

   LD_stats <- list(pair_info=pairs, R2=R2_pairs)
   return(LD_stats)
#plot(as.numeric(pairs[,3]), D_pairs)
#plot(as.numeric(pairs[,3]), R2_pairs)

}


