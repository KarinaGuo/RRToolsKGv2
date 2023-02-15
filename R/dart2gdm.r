#' Prepares input files that are ready
#' for analysis using R package gdm
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop     -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for localdiff analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- dart2gdm(gms, etc)
#' }

dart2gdm <- function(dms, basedir, species, dataset, pop, env_var=FALSE, varfile, climdir) {

   # Step 1, get the genotypes ready
   treatment <- dms$treatment 
   if (dms$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. \n")
   } else {
      cat(" The genotype matrix does not appear to have altcount genotype encoding. \n"); stop()
   }

   # Population allele stats
   population_allele_stats  <- calculate.population.allele.stats(dms, pop)
   population_spatial_dist  <- population.pw.spatial.dist(dms, pop)
   population_allele_fst    <- population.pw.Fst(dms, pop, basedir,species,dataset,maf_val = 0.1, miss_val = 0.1)

   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
      population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
      population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
      population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
      population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
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

   gdm_dir    <- paste(basedir,species,"/popgen/",treatment,"/gdm", sep="")
   
   if(!dir.exists(gdm_dir)) {
      cat("  gdm directory: ", gdm_dir, " does not exist and is being created. \n")
      dir.create(gdm_dir)
   } else {
      cat("  gdm directory: ", gdm_dir, " already exists, content will be overwritten. \n")
   }


   allele_fst_matrix      <- population_allele_fst$Fst

   if ( any(allele_fst_matrix < 0) ) {
      cat("  Warning: Fst value less than zero found, setting to zero for GDM... \n")
      allele_fst_matrix[ allele_fst_matrix < 0 ] <- 0

   }

   long_lat               <- population_spatial_dist$pop_info$lon_lat
   labels                 <- as.matrix(population_spatial_dist$pop_info$names,ncol=1)

   if (env_var) {
      clim                   <- get_climate_data(varfile, climdir, long_lat)
      environ                <- cbind(labels, long_lat, clim)
   } else {
      environ                <- cbind(labels, long_lat)
   }

   colnames(environ)[1]   <- "sites"
   allele_fst_file        <- paste(gdm_dir,"/allele_fst.txt",sep="")
   environ_data_file      <- paste(gdm_dir,"/environ_data.txt",sep="")
 
   write.table(allele_fst_matrix, allele_fst_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")
   write.table(environ, environ_data_file, col.names=TRUE, row.names=FALSE, quote=FALSE,sep=" ")

   return(gdm_dir)

}
