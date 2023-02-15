#' Calculates alternate allele counts and numbers of sampled alleles for  
#' a dart object. Results appropriate for BEDASSLE input 
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param pop -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' popstats <- calculate.population.allele.stats(dart_data, pop)
#' }

calculate.population.allele.stats <- function(dart_data, population) {

   gt      <- dart_data$gt
   num_loc <- ncol(gt)
   meta    <- dart_data$meta
   p       <- population

   pop_info <- get.population.info(meta, p, method="average")

   count  <- mat.or.vec(pop_info$number, num_loc)
   sample <- mat.or.vec(pop_info$number, num_loc)
   freq   <- mat.or.vec(pop_info$number, num_loc)

   for (i in 1:pop_info$number) {
      for (j in 1:num_loc) {
         i_pop_indices  <- which(p == pop_info$names[i])
         g  <- gt[i_pop_indices,j]
         gs <- g[ which(!is.na(g)) ]
         count[i, j]  <- sum(gs)
         sample[i, j] <- length(gs)*2
     
         if (length(gs)*2 > 0) {
            count[i, j]  <- sum(gs)
            sample[i, j] <- length(gs)*2
            freq[i, j]   <- sum(gs) / (length(gs)*2)
         }
         else {
            count[i, j]  <- NA
            sample[i, j] <- NA
            freq[i, j]   <- NA
         }
      }
   }

   minor      <- count
   ind_major  <- which( ( 2*colSums(count) ) > colSums(sample)  )
   minor[ , ind_major ] <- sample[ , ind_major] - count[ , ind_major ]

   colnames(sample) <- colnames(gt)
   colnames(count)  <- colnames(gt)
   colnames(freq)   <- colnames(gt)
   colnames(minor)  <- colnames(gt)
   rownames(sample) <- pop_info$names
   rownames(count) <- pop_info$names
   rownames(freq)   <- pop_info$names
   rownames(minor)  <- pop_info$names

   pop_stats <- list(count=count, sample=sample, minor=minor, freq=freq) 
   return(pop_stats)
}

