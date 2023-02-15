#' Returns a dart data object as a genlight object
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param meta_data -- name of a file containing meta-data for samples  
#' @return a genlight \{adegenet\} object containing the same data
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gl <- dart2gl(dart_data, meta.csv)
#'}

dart2gl <- function(dart_data, basedir, species, dataset, meta_data=FALSE, return_as_file=FALSE) {

   require(adegenet)
   treatment <- dart_data$treatment 

   if (dart_data$encoding == "altcount") {

      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genlight. \n")
      gl_gt <- dart_data$gt

   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   if (!meta_data) {
      cat(" Meta data file not specified \n")
   } else {
      cat(" Ummm...!! \n")
   }

   gl_nuc <- gsub(">", "/", dart_data$locus_nuc)
   
   dart_gl <- new("genlight", gen=gl_gt, ploidy=2, ind.names=dart_data$sample_names, loc.names=dart_data$locus_names, loc.all=gl_nuc)


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

   gl_dir    <- paste(basedir,species,"/popgen/",treatment,"/genlight", sep="")
   
   if(!dir.exists(gl_dir)) {
      cat("  genlight directory: ", gl_dir, " does not exist and is being created. \n")
      dir.create(gl_dir)
   } else {
      cat("  genlight directory: ", gl_dir, " already exists, content will be overwritten. \n")
   }

   gl_file   <- paste(gl_dir,"/",species,"_",dataset,".rda",sep="")

   save(dart_gl, file=gl_file)

   if (return_as_file) { 
      return(gl_file)
   } else {
      return(dart_gl)
   }
}
