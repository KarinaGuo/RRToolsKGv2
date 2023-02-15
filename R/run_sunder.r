#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param sun_par   -- list of sunder params  [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param name      -- arbitrary name for the dataset [required]
#' @param treatment -- the treatment string, which is used to find file [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_sunder(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 

run_sunder <- function(s, basedir, species, dataset, treatment) {

   require(Sunder)
   s_MCMCCV  <- MCMCCV(gen=s$gen, D_G=s$D_G, D_E=s$D_E, nit=s$nit, thinning=s$thinning, theta.max=s$theta.max, theta.init=s$theta.init, run=c(FALSE,TRUE,FALSE), ud=s$ud, n.validation.set=s$n.validation.set, print.pct=TRUE)
   

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

   su_dir    <- paste(basedir,species,"/popgen/",treatment,"/sunder", sep="")
   
   if(!dir.exists(su_dir)) {
      cat("  sunder directory: ", su_dir, " does not exist and is being created. \n")
      dir.create(su_dir)
   } else {
      cat("  sunder directory: ", su_dir, " already exists, content will be overwritten. \n")
   }

   su_run_file   <- paste(su_dir,"/",species,"_",dataset,"_RUN_udall.rda",sep="")

   save(s_MCMCCV,file=su_run_file)

}



