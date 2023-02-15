#' Returns a matrix of spatial distances among populations
#' 
#' Input is a dart meta data object
#' 
#' @param dart_data  -- a dart meta data object [Required]
#' @param population -- a dart meta data object [Required]
#' @return Fst estimates across a nominated set of populations
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

multipop.Fst <- function(dart_data, population, basedir, species, dataset, maf_val=0.2, miss_val=0.2) {

   meta <- dart_data$meta
   p    <- population
   pop_info <- get.population.info(meta, p, method="average")

   require(SNPRelate)
   gds_file <- dart2gds(dart_data, basedir, species, dataset)
   gds      <- snpgdsOpen(gds_file)

   fst      <- snpgdsFst(gds, population=as.factor(p), method="W&H02", sample.id=dart_data$sample_names, maf=maf_val, missing.rate=miss_val, with.id=TRUE)

   snpgdsClose(gds)
   return(fst)
}

