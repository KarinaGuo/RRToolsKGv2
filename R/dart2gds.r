#' Returns a dart data object as a genind object
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param meta_data -- name of a file containing meta-data for samples  
#' @return a genind \{adegenet\} object containing the same data
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gds <- dart2gds(dart_data, basedir, species, dataset, meta_data=FALSE)
#' }

dart2gds <- function(dart_data, basedir, species, dataset, meta_data=FALSE) {

   treatment <- dart_data$treatment 

   if (dart_data$encoding == "altcount") {

      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to gds. \n")
      
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   if (!meta_data) {
      cat(" Meta data file not specified \n")
   } else {
      cat(" Ummm...!! \n")
   }


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

   gds_dir    <- paste(basedir,species,"/popgen/",treatment,"/gds", sep="")
   
   if(!dir.exists(gds_dir)) {
      cat("  gds directory: ", gds_dir, " does not exist and is being created. \n")
      dir.create(gds_dir)
   } else {
      cat("  gds directory: ", gds_dir, " already exists, content will be overwritten. \n")
   }

   gds_file   <- paste(gds_dir,"/",species,"_",dataset,".gds",sep="")

   snpgdsCreateGeno(gds_file, genmat = dart_data$gt, sample.id = dart_data$sample_names, snp.id = colnames(dart_data$gt), snpfirstdim=FALSE)

   return(gds_file)

}
