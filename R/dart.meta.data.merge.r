#' Merges dart and meta data objects into a single object 
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param meta_data -- a meta data object [Required]  
#' @return a dart_meta_data object
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

dart.meta.data.merge <- function(dart_data, meta_data) {

   d <- dart_data
   m <- meta_data

   # find overlapping samples in dart and meta datasets
   dart_meta_int <- intersect(rownames(d$gt), m$sample_names)
   id <- which(rownames(d$gt) %in% dart_meta_int)
   im <- which(m$sample_names %in% dart_meta_int)

   d$gt             <- d$gt[ id , ]
   d$sample_names   <- d$sample_names[ id ]

   m$analyses       <- m$analyses[ im, , drop=FALSE]
   m$lat            <- m$lat[ im ]
   m$long           <- m$long[ im ]
   m$sample_names   <- m$sample_names[ im ]
   m$site           <- m$site[ im ]

   # reorder the dart (genotype) data to be consistent 
   # consistent with sample order in meta data

   oid <- order(match(rownames(d$gt),m$sample_names))
   d$gt             <- d$gt[ oid , ]
   d$sample_names   <- d$sample_names[ oid ]

   # remove any fixed snps
   dp <- remove.fixed.snps(d)

   # check sample names are the same
   if ( !identical(dp$sample_names, m$sample_names) ) {
      cat(" Fatal error: after match and order, samples not identical in dart and meta \n"); stop()
   } else {
      cat(" Samples successfully matched and ordered: merging dart and meta objects \n")
   }


   dart_data_meta <- dp
   dart_data_meta$meta <- m

   return(dart_data_meta)

}
