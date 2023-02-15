#' Recodes genotypes supplied from DArT format to count of alternate allele
#' Also changes the status from the $encoding element of dart data list 
#' from 'DArT' to 'altcount'

#'
#' @param dart_data -- dart data list  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' encode.dart2altcount(dart_data)


encode.dart2altcount <- function(dart_data) {

     dart_data_proc <- dart_data
     
     if (dart_data$encoding == "DArT") {
          cat(" DArT object found with dart encoding (0=hom ref; 1=hom alt; 2=het) \n");
          cat(" Rewriting with altcount encoding     (0=hom ref; 1=het; 2=hom alt) \n");  

          dart_data_proc$encoding <- "altcount"

          dart_data_proc$gt[  which( dart_data$gt == 1) ] <- 2
          dart_data_proc$gt[  which( dart_data$gt == 2) ] <- 1

          return(dart_data_proc)

     } else {
          cat(" Fatal Error: data not in dart format \n"); stop()
     }

}
