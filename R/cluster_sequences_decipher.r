#' performs sequence clustering using package decipher
#' 
#' @return cdhit output
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' ld_files <- dart2localdiff(gms, etc)
#' }


   cluster_sequences_decipher <- function(fa_fil, cl_fil, cutoff=cutoff, method="complete",np=NULL) {

      library(DECIPHER)
      fas <- fa_fil
      dna <- readDNAStringSet(fas)

      if (method == "inexact") {
         idc <- IdClusters(myXStringSet=dna, method="inexact", cutoff=cutoff)
      }

      if (method == "complete") {
         d <- DistanceMatrix(dna, processors=np) # returns an object of class 'dist'
         idc <- IdClusters(d, method="complete", cutoff=cutoff, processors=np)
      }

      write.table(idc,cl_fil,col.names=FALSE, row.names=TRUE,quote=FALSE,sep=",")

   }
