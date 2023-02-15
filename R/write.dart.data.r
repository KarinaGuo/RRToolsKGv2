#' Write a dart data list object to an R workspace file  
#'
#' write.dart.data() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param dart_data -- dart data list  [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param name      -- arbitrary name for the dataset [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' write.dart.data(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 
write.dart.data <- function(dart_data, basedir, species, dataset, meta=FALSE) {

   name <- dataset
   if (meta) {
      treatment <- dart_data$dart_data$treatment 
   } else {
      treatment <- dart_data$treatment 
   }

   

   dir <- paste(basedir, species, "/dart_standard",sep="")

   if(!dir.exists(dir)) {
      cat("  Dart standard directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Dart standard directory: ", dir, " already exists... content might be overwritten. \n")
   }

   ds_dir <- paste(basedir, species, "/dart_standard/",treatment,sep="")

   if(!dir.exists(ds_dir)) {
      cat("  Standard directory: ", ds_dir, " does not exist and is being created. \n")
      dir.create(ds_dir)
   } else {
      cat("  Standard directory: ", ds_dir, " already exists, content will be overwritten. \n")
   }

   ds_file   <- paste(ds_dir,"/",species,"_",dataset,".rda",sep="")

   save(dart_data, file=ds_file)

   return(ds_file)

}



