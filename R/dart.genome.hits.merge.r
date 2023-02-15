#' Merges dart and meta data objects into a single object 
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param basedir -- directory for analysis [Required]  
#' @param species -- the 8 letter species identifier [Required]  
#' @param dataset -- name of dataset [Required]  
#' @return dart_data object with genome info fields
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_data <- dart.genome.hits.merge(dart_data, basedir, species, dataset)
#'}

dart.genome.hits.merge <- function(dart_data, basedir, species, dataset) {

   d <- dart_data

   genome_file <- paste(basedir, species, "/genome/",species,"_",dataset,"_aln_eucagen.rda",sep="")

   if(!file.exists(genome_file)) {
      cat(" Fatal error: File containing genomic locations does not exist."); stop();
   } 
 
   load(genome_file)

   g <- eucagenhits

   # find overlapping samples in dart and meta datasets
   dart_genome_int <- intersect(colnames(d$gt), g$CloneID)
   id <- which(colnames(d$gt) %in% dart_genome_int)
   ig <- which(g$CloneID %in% dart_genome_int)

   g  <- g[ ig, ]

   # reorder the genome location data to be consistent 
   # with locus order in the genotype 

   oid <- order(match(g$CloneID, rownames(d$gt)))
   g   <- g[ oid , ]

   # check sample names are the same
   if ( !identical(colnames(d$gt), g$CloneID) ) {
      cat(" Fatal error: after match and order, samples not identical in dart and genome information \n"); stop()
   } else {
      cat(" Samples successfully matched and ordered: merging dart object and genome info \n")
   }

   dart_data_proc <- d
   dart_data_proc$eucagen_locations <- g

   return(dart_data_proc)

}
