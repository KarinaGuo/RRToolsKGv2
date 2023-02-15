#' Returns a list of allele frequencies and heterozygosities 
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @return a list containing a vector of alternate allele frequencies
#' and a vector of observed heterozygosities
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' AF.summary <- calculate.AF.summary(dart_data)


calculate.AF.summary <- function(dart_data) {

   treatment <- dart_data$treatment 

   if (dart_data$encoding == "altcount") {

      # cat(" Dart data object for ", dataset, "in species", species, "\n")

   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

 
   x <- dart_data$gt

   f0hom <- function(x) {
      length(which(x == 0))
   }

   fhet <- function(x) {
      length(which(x == 1))
   }

   f1hom <- function(x) {
      length(which(x == 2))
   }

   hom0 <- apply(x, 2, f0hom)
   het  <- apply(x, 2, fhet)
   hom1 <- apply(x, 2, f1hom)
 
   P <- (2*hom1 + het) / (2 * (hom0 + het + hom1))
   H <- het / (hom0 + het + hom1)

   AF.summary <- list(P=P, H=H)

   return(AF.summary)

}
