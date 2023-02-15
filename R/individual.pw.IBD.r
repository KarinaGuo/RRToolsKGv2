#' Returns a matrix of spatial distances among populations
#' 
#' Input is a dart meta data object
#' 
#' @param dart_data -- a dart meta data object [Required]
#' @return a matrix containing spatial distances between populations
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

individual.pw.IBD <- function(dart_data, basedir, species, dataset, maf=0.1) {

   meta <- dart_data$meta

   require(SNPRelate)
   gds_file <- dart2gds(dart_data, basedir, species, dataset)
   gds <- snpgdsOpen(gds_file)

   IBD <- snpgdsIBDMoM(gds, maf=maf, missing.rate=0.2, num.thread=1, kinship=TRUE)

   snpgdsClose(gds)
   return(IBD)
}
