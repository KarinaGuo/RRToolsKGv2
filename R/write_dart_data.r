#' Write a dart data list object to an R workspace file  
#'
#' @param dart_data -- dart data list  [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param dataset   -- arbitrary name for the dataset [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' write.dart.data(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 
write_dart_data <- function(dart_data, basedir, species, dataset) {
 
   if(!dir.exists(paste(basedir, species, sep=""))) {
      cat("  Analysis directory for this species does not exist and is being created at ", paste(basedir, species, sep=""), "\n" )
      dir.create(paste(basedir, species, sep=""))
   } else {
      
   }



   dir <- paste(basedir, species, "/dart_standard",sep="")

   if(!dir.exists(dir)) {
      cat("  Dart standard directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Dart standard directory: ", dir, " already exists... content might be overwritten. \n")
   }

   treatment <- dart_data$treatment
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



