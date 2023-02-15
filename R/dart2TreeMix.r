#' Returns a dart data object in an object that is ready
#' for analysis using software SpaceMix
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop     -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2rangexp(dms, meta.csv)
#' }

dart2TreeMix <- function(dms, basedir, species, dataset, pop) {

   # Step 1, get the genotypes ready
   treatment <- dms$treatment 
   if (dms$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # Population allele stats
   population_allele_stats  <- calculate.population.allele.stats(dms, pop)
   
   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
      population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
      population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
      population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
      population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
   }

   # make directory, write files 
   dir <- paste(basedir, species, "/popgen",sep="")
   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
   }

   dir <- paste(basedir, species, "/popgen/",treatment,sep="")

   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists...  \n")
   }

   tm_dir    <- paste(basedir,species,"/popgen/",treatment,"/TreeMix", sep="")
   
   if(!dir.exists(tm_dir)) {
      cat("  treemix directory: ", tm_dir, " does not exist and is being created. \n")
      dir.create(tm_dir)
   } else {
      cat("  treemix directory: ", tm_dir, " already exists, content will be overwritten. \n")
   }

   tm_object_file   <- paste(tm_dir,"/",species,"_",dataset,"tm_in",sep="")

   counts       <- population_allele_stats$minor
   sample_sizes <- population_allele_stats$sample
   counts_alt   <- sample_sizes - counts

   ref_alt    <- paste(counts, counts_alt, sep=",")
   tree_mix   <- t(matrix(ref_alt,nrow=nrow(counts)))
   colnames(tree_mix) <- rownames(population_allele_stats$sample)
   
   write.table(tree_mix, tm_object_file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")

   return(tm_object_file)

}
