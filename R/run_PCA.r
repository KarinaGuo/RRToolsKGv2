#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param gl_file   -- location of an R file containing genlight object  [required]
#' @param basedir   -- Base dirctory  [required]
#' @param species   -- species name  [required]
#' @param dataset   -- dataset identifier  [required]
#' @param treatment -- treatment  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_PCA(gl_file, basedir, species, dataset, treatment)
#
 

run_PCA <- function(gl_file, basedir, species, dataset, treatment) {

   require(adegenet)

   if (file.exists(gl_file)) {

      cat("opening ", gl_file, "for PCA \n")
      load(gl_file)   

   } else {
      cat("Fatal error: file with genlight object ", gl_file, " does not exist \n"); stop();
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

   gl_dir    <- paste(basedir,species,"/popgen/",treatment,"/genlight", sep="")
   
   if(!dir.exists(gl_dir)) {
      cat("  genlight directory: ", gl_dir, " does not exist and is being created. \n")
      dir.create(gl_dir)
   } else {
      cat("  genlight directory: ", gl_dir, " already exists, content will be overwritten. \n")
   }

   pca_file   <- paste(gl_dir,"/",species,"_",dataset,"_PCA.rda",sep="")
   pca_plot   <- paste(gl_dir,"/",species,"_",dataset,"_PCA.pdf",sep="")

   gl_pca <- glPca(dart_gl, nf=5, parallel=FALSE)

   pdf(file=pca_plot)
      scatter(gl_pca)
   dev.off()

   save(gl_pca, file=pca_file)
}



