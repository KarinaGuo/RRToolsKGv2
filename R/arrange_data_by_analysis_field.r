#' Merges dart and meta data objects into a single object 
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a merged dart-meta data object [Required]
#' @param field     -- vector data field names for splitting [Required]  
#' @param basedir   -- name of the base directory [Required]  
#' @param species   -- species name [Required]  
#' @param dataset   -- dataset identifier [Required]
#' @return object   -- a list of file names for the data objects that are produced
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

arrange_data_by_analysis_field <- function(dart_data, field, basedir, species, dataset) {

   dm <- dart_data

   for (i in 1:length(dm$meta$analyses)) { 
      dm$meta$analyses[[i]][ dm$meta$analyses[[i]] == ""  ]  <- NA 
      dm$meta$analyses[[i]][ dm$meta$analyses[[i]] == "-"  ] <- NA
   }

   ind_remove      <- which( is.na( dm$meta$analyses[[field]] ) )   
   nm              <- list()

   if (length(ind_remove) > 0) {
         if (is.list(dm$meta$analyses)) {
               nm$analyses     <- lapply(dm$meta$analyses, function(x) x <- x[ -ind_remove]) 
         }
         
         
         nm$lat          <- dm$meta$lat[ -ind_remove ]
         nm$long         <- dm$meta$long[ -ind_remove ]
         nm$sample_names <- dm$meta$sample_names[ -ind_remove ]
         nm$site         <- dm$meta$site[ -ind_remove ]
         ndm             <- merge_gt_meta_data(dm,nm)

   } else {
         nm$analyses     <- dm$meta$analyses 
         nm$lat          <- dm$meta$lat
         nm$long         <- dm$meta$long
         nm$sample_names <- dm$meta$sample_names
         nm$site         <- dm$meta$site
         ndm             <- merge_gt_meta_data(dm,nm)
   }

   ndm$treatment <- paste(ndm$treatment,"_F", field, sep="")
   file          <- write_dart_data(ndm, basedir, species, dataset)

   return(ndm)

}


