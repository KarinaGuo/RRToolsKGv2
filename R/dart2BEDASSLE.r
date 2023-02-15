#' Returns a dart data object in an object that is ready
#' for analysis using software BEDASSLE
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

dart2BEDASSLE <- function(dms, basedir, species, dataset, pop) {

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
   population_spatial_dist  <- population.pw.spatial.dist(dms, pop)

   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
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

   bd_dir    <- paste(basedir,species,"/popgen/",treatment,"/bedassle", sep="")
   
   if(!dir.exists(bd_dir)) {
      cat("  bedassle directory: ", bd_dir, " does not exist and is being created. \n")
      dir.create(bd_dir)
   } else {
      cat("  bdeassle directory: ", bd_dir, " already exists, content will be overwritten. \n")
   }

   bd_object_file   <- paste(bd_dir,"/",species,"_",dataset,".rda",sep="")

   counts       <- population_allele_stats$count
   sample_sizes <- population_allele_stats$sample
   D            <- population_spatial_dist$S/100000
   E            <- population_spatial_dist$S
   E[,]         <- 0
   #E[1:2,1:2]   <- 0.1
   k            <- nrow(counts)
   loci         <- ncol(counts)
   delta        <- 0.001
   aD_stp       <- 0.002
   aE_stp       <- 0.001
   a2_stp       <- 0.002
   thetas_stp   <- 0.05
   mu_stp       <- 0.17
   ngen         <- 10000000
   printfreq    <- ngen/200
   savefreq     <- ngen/100
   samplefreq   <- 100
   directory    <- bd_dir
   prefix       <- paste(species,"_",dataset,sep="")

   bd <- list( counts=counts, sample_sizes=sample_sizes, D=D, E=E, k=k, loci=loci, delta=delta, aD_stp=aD_stp, aE_stp=aE_stp, a2_stp=a2_stp, thetas_stp=thetas_stp, mu_stp=mu_stp, ngen=ngen, printfreq=printfreq, savefreq=savefreq, samplefreq=samplefreq, directory=directory, prefix=prefix)
   save(bd, file=bd_object_file)

   return(bd)

}
