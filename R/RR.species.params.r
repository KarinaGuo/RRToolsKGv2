#' Read from a table of RR genotype datasets, get species-specific parameters
#'
#' @param paramdir     -- name of directory containing a file with a standard parameter set [required]
#' @param species_file -- name of the file containing information about the genotype data [required]
#' @return             -- A list of parameters, which can be used in various data filtering functions
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' \dontrun{
#' p <- RR.species.params(RandRparam, species_file)
#' } 

RR.species.params <- function(paramdir, species_file, target_species="all") {

   file <- paste(paramdir, species_file, sep="")

   cat("\n")
   cat(" Opening the list of species and parameters:", file,"\n")

   species_info <- read.table(file, sep="\t", header=TRUE)


   if (target_species == "all") {
      cat(" Reading information for all species that are listed\n")
   } else {
      ind_species <- which(species_info$species %in% target_species)
      species_info <- species_info[ ind_species, ]
   }

   num_species  <- length(species_info[,"species"])
   s <- list()

   for (i in 1:num_species) {
      species  <- as.character(species_info[i,"species"])
      topskip  <- as.numeric(species_info[i,"topskip"])
      nmetavar <- as.numeric(species_info[i,"nmetavar"])
      dataset  <- as.character(species_info[i,"dataset"])
      genus_spp  <- as.character(species_info[i,"genus_spp"])
      tmp <- list(species=species, topskip=topskip, nmetavar=nmetavar, dataset=dataset, genus_spp=genus_spp)
      s[[species]] <- tmp
   }
   return(s)
}
