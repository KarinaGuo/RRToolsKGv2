#' Returns a dart data object in appropriate files for analysis 
#' with the Peter and Slatkin method of detecting a range expansion 
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2rangexp(dart_data, meta.csv)
#' }

dart2rangexp <- function(dart_data, basedir, species, dataset, pop) {

   # Step 1, get the genotypes ready
   treatment <- dart_data$treatment 
   if (dart_data$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # genotypes -- change char for missing
   gt <- dart_data$gt
   gt[ is.na(gt) ] <- "?"

   # meta -- make a table
   mt <- cbind(dart_data$meta$sample_names, dart_data$meta$lat, dart_data$meta$long, pop)
   colnames(mt) <- c("id","latitude","longitude","pop")

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

   re_dir    <- paste(basedir,species,"/popgen/",treatment,"/rangexp", sep="")
   
   if(!dir.exists(re_dir)) {
      cat("  RE directory: ", re_dir, " does not exist and is being created. \n")
      dir.create(re_dir)
   } else {
      cat("  RE directory: ", re_dir, " already exists, content will be overwritten. \n")
   }

   re_genotype_file   <- paste(re_dir,"/",species,"_",dataset,".snapp",sep="")
   re_spatial_file    <- paste(re_dir,"/",species,"_",dataset,".space",sep="")

   write.table(gt, file=re_genotype_file, quote=FALSE, row.names=TRUE, col.names=FALSE, sep=",")
   write.table(mt, file=re_spatial_file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")

   re_files <- list(re_genotype_file=re_genotype_file, re_spatial_file=re_spatial_file)

   return(re_files)

}
