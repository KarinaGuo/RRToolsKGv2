#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param lea_file   -- location of lea file  [required]
#' @param basedir   -- Base dirctory  [required]
#' @param species   -- species name  [required]
#' @param dataset   -- dataset identifier  [required]
#' @param treatment -- treatment  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_PCA(lea_file, basedir, species, dataset, treatment)
#
 

run_SNMF <- function(lea_file, basedir, species, dataset, treatment) {

   require(LEA)

   if (file.exists(lea_file)) {

      cat("LEA file", lea_file, "exists, preparing to run \n")
 
   } else {
      cat("Fatal error: the LEA file ", lea_file, " does not exist \n"); stop();
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

   out_file   <- paste(lea_dir,"/",species,"_",dataset,"_LEA.rda",sep="")
   lea_plot   <- paste(lea_dir,"/",species,"_",dataset,"_LEA_K.pdf",sep="")

   snmfK20R10=snmf(lea_file, K=1:20, entropy = TRUE, repetitions = 10, project = "new")
   pdf(file=lea_plot)
   plot(snmfK20R10, lwd = 5, col = "red", pch=1)
   dev.off()

   save(snmfK20R10, file=out_file)

   return(out_file)
}



