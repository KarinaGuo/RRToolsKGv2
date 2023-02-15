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
#' dart_struct <- dart2struct(dart_data, pop_prior=FALSE)
#' }

dart2struct <- function(dart_data, basedir, species, dataset, use_pops=FALSE, pops) {

   if (dart_data$encoding == "altcount") {

      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to STRUCTURE format. \n")
      st_gt <- dart_data$gt

   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }



   st_gt[ is.na(st_gt) ] <- "-9 -9"
   st_gt[ which(st_gt == 0) ] <- "0 0"
   st_gt[ which(st_gt == 1) ] <- "0 1"
   st_gt[ which(st_gt == 2) ] <- "1 1"

   treatment <- dart_data$treatment 

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

   str_dir    <- paste(RandRbase,species,"/popgen/",treatment,"/struct", sep="")
   
   if(!dir.exists(str_dir)) {
      cat("  structure directory: ", str_dir, " does not exist and is being created. \n")
      dir.create(str_dir)
   } else {
      cat("  structure directory: ", str_dir, " already exists, content will be overwritten. \n")
   }

   if (!use_pops) {
      cat(" Population prior info is not specified... \n")
   } else {
      cat(" Adding population info to data matrix \n")
      rownames(st_gt) <- paste(rownames(st_gt), pops, sep=" ")

   }

   str_file   <- paste(str_dir,"/",species,"_",dataset,".struct",sep="")
   write.table(st_gt, file=str_file,sep=" ",quote=FALSE, row.names = TRUE, col.names = FALSE)

   return(str_file)

}
