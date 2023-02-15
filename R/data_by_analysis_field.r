#' Merges dart and meta data objects into a single object 
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a merged dart-meta data object [Required]
#' @param fields    -- vector data field names for splitting [Required]  
#' @param basedir   -- name of the base directory [Required]  
#' @param species   -- species name [Required]  
#' @param dataset   -- dataset identifier [Required]
#' @param object    -- if FALSE, return a vector of file names    
#' @return object   -- a list of file names for the data objects that are produced
#' @return field    -- the data object for the nominated field
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

data.by.meta.fields <- function(dart_data, fields, basedir, species, dataset, object="files") {

   dm <- dart_data
   f  <- fields

   num_fields <- length(f)

   file_vector <- NULL

   for (i in 1:num_fields) {

      field <- fields[i]

      dm$meta$analyses[dm$meta$analyses == ""]  <- NA
      dm$meta$analyses[dm$meta$analyses == "-"] <- NA

      ind_remove      <- which( is.na( dm$meta$analyses[ , field ] ) )   
      nm              <- list()

      if (length(ind_remove) > 0) {
         nm$analyses     <- dm$meta$analyses[ -ind_remove,  , drop=FALSE] 
         nm$lat          <- dm$meta$lat[ -ind_remove ]
         nm$long         <- dm$meta$long[ -ind_remove ]
         nm$sample_names <- dm$meta$sample_names[ -ind_remove ]
         nm$site         <- dm$meta$site[ -ind_remove ]
         ndm             <- dart.meta.data.merge(dm,nm)

      } else {
         nm$analyses     <- dm$meta$analyses 
         nm$lat          <- dm$meta$lat
         nm$long         <- dm$meta$long
         nm$sample_names <- dm$meta$sample_names
         nm$site         <- dm$meta$site
         ndm             <- dart.meta.data.merge(dm,nm)
      }

      ndm$treatment <- paste(ndm$treatment,"_Field_", field, sep="")
      file <- write.dart.data(ndm, basedir, species, dataset)
      file_vector <- c(file_vector, file)

      if (object==field) { dm_return <- ndm }

   }

   if (object=="files") { return(file_vector) }
   else { 
      if (exists("dm_return")) { return(dm_return) }
      else {
         cat("Fatal error: nominate to return ", object, "but not found\n"); stop()
      } 
   }
}


