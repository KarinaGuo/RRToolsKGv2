#' Returns a matrix of spatial distances among individuals
#' 
#' Input is a dart meta data object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a matrix containing spatial distances between populations
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

individual.pw.spatial.dist <- function(dart_data) {

   require(geosphere)

   meta <- dart_data$meta

   sample_number <- length(meta$sample_names)

   S <- mat.or.vec(sample_number, sample_number)

   for (i in 1:sample_number) {
      for (j in 1:sample_number) {

         if (i > j) {
            LLi <- c(meta$long[i], meta$lat[i])
            LLj <- c(meta$long[j], meta$lat[j])
            Dij <- distCosine(as.numeric(LLi),as.numeric(LLj))
            S[i,j] <- Dij
            S[j,i] <- Dij
         }
      }
   }

   colnames(S) <- meta$sample_names
   rownames(S) <- meta$sample_names
   slist <- list(S=S) 
   return(slist)
}
