#' Import metadata for a DArT dataset 
#'
#' read.meta.data() opens a meta data file associated with a dart dataset 
#'
#' @param datafile -- name of xls file containing the meta data [required]
#' @param nas      -- missing data character [default "-"]
#' @return         -- a dart meta data object, consisting of sample names
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' meta_data <- read.dart.meta(datafile=, nas="-")
#
 
read.meta.data <- function(dart_data, basedir, species, dataset, fields=1, nas="-") {

   metafile <- paste(basedir, species, "/meta/", species, "_", dataset,"_meta.xlsx", sep="")
   if ( file.exists(metafile) ) {
      cat("\n")
      cat(" Reading data file:", metafile,"\n")
   } else {
      cat(" Fatal Error: the metadata file", metafile,"does not exist \n"); stop()
   }

   ### Open and parse meta data 
   m  <- read_excel(metafile, sheet=1, col_names=TRUE, na=nas)

   if ( any( names(m) == "sample" ) & any(names(m) == "lat") & any(names(m) == "long") ) {
      cat(" Found sample, lat and long columns in metadata \n")
   } else {
      cat(" Fatal Error: did not find important sample, lat or long column in metadata \n"); stop()
   }

   mm  <- match(rownames(dart_data$gt), m$sample)
   mi  <- intersect(m$sample, rownames(dart_data$gt))
   #mm_real <- mm[which(!is.na(mm))]
   mm_real <- which(m$sample %in% mi)

   num_dart_samples  <- nrow(dart_data$gt) 
   num_meta_samples  <- nrow(m)
   num_match_samples <- length( mm_real )


   if (num_match_samples > 1) {

      cat(" Found metadata for ", num_meta_samples, " samples \n")
      cat(" This includes overlap with ", num_match_samples, " samples \n")
      cat(" out of ", num_dart_samples, "in DArT genotypes \n")

   } else {
      cat(" Fatal Error: no matching sample information between meta-data and DArT genotypes \n"); stop()
   }

   missing_in_meta <- rownames(dart_data$gt)[is.na(m$sample[mm])]

   meta_ordered <- m[ mm_real, ]
   sample_names <- as.character(meta_ordered$sample)
   site         <- as.character(meta_ordered$site)
   lat          <- as.numeric(meta_ordered$lat)
   long         <- as.numeric(meta_ordered$long)

   analyses <- "none"

   if ( fields == 1 ) {
      cat(" Adding analysis field to meta data list \n") 
      analyses <- matrix(meta_ordered[, 5],)
      colnames(analyses) <- colnames(meta_ordered)[5]
   }


   if ( fields > 1 ) {
      cat(" Adding analysis fields to meta data list \n") 
      analyses <- as.matrix(meta_ordered[, 5:(5+fields-1)] )

   }

   meta_data <- list(sample_names=sample_names, site=site, lat=lat, long=long, analyses=analyses)

   return(meta_data)
}


