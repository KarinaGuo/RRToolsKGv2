#' Recodes genotypes supplied from count of alternate allele to count of minor allele
#' Also changes the status from the $encoding element of data list 

#'
#' @param dart_data -- dart data list  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' encode.dart2altcount(dart_data)


encode.altcount2mincount <- function(dart_data) {

     dart_data_proc <- dart_data
     
     if (dart_data$encoding == "altcount") {
          cat(" DArT object found with altcount encoding (0=hom ref; 1=het; 2=hom alt) \n");
          cat(" Rewriting with mincount encoding     (0=hom maj; 1=het; 2=hom min) \n");  

          dart_data_proc$encoding <- "mincount"

          afsumm <- calculate.AF.summary(dart_data)
          indmaj <- which(afsumm$P > 0.5)

          for (i in indmaj) {
             dart_data_proc$gt[  which( dart_data$gt[,i] == 0),i ] <- 2
             dart_data_proc$gt[  which( dart_data$gt[,i] == 2),i ] <- 0
          } 
 
          return(dart_data_proc)

     } else {
          cat(" Fatal Error: data not in dart format \n"); stop()
     }

}
