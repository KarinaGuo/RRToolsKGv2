#' Prepares SNMF Q matrices to act as explanator variables in GDM
#'
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param treatment  
#' @param process [TRUE / FALSE] -- indicates whether to process Q matrix, 
#'                                  or if it has been processed previously 
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- SNMF2gdm(basedir, species, dataset, treatment, process=TRUE, Ksel=2)
#' }

SNMF2gdm <- function(basedir, species, dataset, treatment, process=FALSE, pop=NULL, Ksel=NULL) {

   if (process) {
      if( is.null(pop) | is.null(Ksel) ) {
         cat( "   Attempting to process SNMF output, but pop info or Ksel not supplied ...\n"); stop();
      }
      pop_Q_file <- process_SNMF_for_gdm(basedir, species, dataset, treatment, pop, Ksel)

   } else {
      pop_Q_file  <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,"_popQ.txt", sep="")
   }

   gdm_E_file  <- paste(basedir,species,"/popgen/",treatment,"/gdm/environ_data.txt", sep="")
   gdm_EQ_file <- paste(basedir,species,"/popgen/",treatment,"/gdm/environ_Q_data.txt", sep="")

   if (file.exists(pop_Q_file)) { 
      cat(   "   Found population Q values ...\n") 
   } else { 
      cat(   "   Cannot find population Q files.  \n"); stop() 
   }
   
   if (file.exists(gdm_E_file)) { 
      cat(   "   Found population E values ...\n") 
   } else { 
      cat(   "   Cannot find gdm environmental variable file. Exiting. \n"); stop() 
   }

   environ  <- read.table(gdm_E_file, header=TRUE, sep=" ") 
   ancestr  <- read.table(pop_Q_file, header=TRUE, sep=" ") 

   index    <- 1:nrow(environ) 
   environ  <- cbind(environ, index)
   environQ <- merge(environ, ancestr, by="sites")
   environQ <- environQ[ order(environQ[,"index"]),]
   environQ <- environQ[,-which( colnames( environQ ) =="index")]

   write.table(environQ, gdm_EQ_file, quote=FALSE, col.names=TRUE, row.names=FALSE)

   return(gdm_EQ_file)

}

