#' parses cdhit output in R
#' 
#' @return cdhit output
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' ld_files <- dart2localdiff(gms, etc)
#' }

parse_cdhit_file <- function(cdh_out) {

   co      <- readLines(cdh_out)
   locus   <- NULL
   cluster <- NULL
   cn  <- 0 
  for (i in 1:length(co)) {

      lin <- co[i]
      
      if (substr(lin,1,1) == ">") { cn <- cn + 1 }
      if (substr(lin,1,1) != ">") { 

         tmpbit  <- strsplit(co[i],">")[[1]][2]
         loct    <- strsplit(tmpbit,"[...]")[[1]][1]
         cluster <- c(cluster, cn)
         locus   <- c(locus, loct)
         
      }
   }

   cdh.df <- data.frame(locus,cluster)
   return(cdh.df)

}
