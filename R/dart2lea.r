#' Writes a dart data object to an lfmm genotype file
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_lfmm <- dart2lfmm(dart_data, meta=FALSE)
#' }

dart2lea <- function(dart_data, basedir, species, dataset, meta_data=FALSE) {

   if (dart_data$encoding == "altcount") {

      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to lfmm. \n")
      lf_gt <- dart_data$gt

   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   if (!meta_data) {
      cat(" Meta data file not specified \n")
   } else {
      cat(" Ummm...!! \n")
   }

   lf_gt[ is.na(lf_gt) ] <- 9

   treatment <- dart_data$treatment 

   dir <- paste(basedir, species,sep="")
   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
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

   lea_dir    <- paste(basedir,species,"/popgen/",treatment,"/lea", sep="")
   
   if(!dir.exists(lea_dir)) {
      cat("  LEA directory: ", lea_dir, " does not exist and is being created. \n")
      dir.create(lea_dir)
   } else {
      cat("  LEA directory: ", lea_dir, " already exists, content will be overwritten. \n")
   }

   lea_file   <- paste(lea_dir,"/",species,"_",dataset,".lfmm",sep="")
   write.table(lf_gt, file=lea_file,sep=" ",quote=FALSE, row.names = FALSE, col.names = FALSE)

   return(lea_file)

}
