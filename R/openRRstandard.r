#' Write a dart data list object to an R workspace file  
#'
#' write.dart.data() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param dart_data -- dart data list  [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param name      -- arbitrary name for the dataset [required]
#' @param treatment -- the treatment string, which is used to find file [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' write.dart.data(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 
openRRstandard <- function(basedir, species, dataset, treatment) {

   ds_dir <- paste(basedir, species, "/dart_standard/",treatment,sep="")
 
   ds_file   <- paste(ds_dir,"/",species,"_",dataset,".rda",sep="")

  if(!file.exists(ds_file)) {
      cat("  Standard data file: ", ds_file, " does not exist  \n"); stop()
   } else {
      cat("  Standard data file: ", ds_file, " will be opened. \n")
   }

   load(ds_file)

   return(dart_data)

}



