#' Processes SNMF output for use as explanatory variables in GDM
#'
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param treatment  
#' @param pop     -- a vector of population assignments
#' @param Ksel    -- the selected number of ancestral populations  
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- SNMF2gdm(basedir, species, dataset, treatment, process=TRUE, Ksel=2)
#' }

process_SNMF_for_gdm <- function(basedir, species, dataset, treatment, pop, Ksel) {

   require(LEA)
   snmf_file    <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,".snmfProject", sep="")
   if (file.exists(snmf_file)) { 
      cat(   "   Calculating population Q values ...\n") 
   } else { 
      cat(   "   SNMF file ", snmf_file," not found. Exiting. \n"); stop() 
   }
 
   snmf_project <- load.snmfProject(snmf_file)
   ce           <- cross.entropy(snmf_project, K = Ksel)
   Rbest        <- which.min(ce)

   cat(   "   Proceeding with sNMF K ", Ksel, "and  run ", Rbest ,"\n")

   qmatrix = Q(snmf_project, K = Ksel, run=Rbest)

   pop_Q_vals <- mat.or.vec(length(unique(pop)), Ksel)

   for (i in 1:Ksel) {
          
          tmp_pop_Q <- sapply(split(qmatrix[,i], pop), mean)
          pop_Q_vals[,i] <- tmp_pop_Q
          if (i==1) { row.names(pop_Q_vals) <- names(tmp_pop_Q) }
   }
    
   transform_pop_Q <- function(pop_Q_vals) {

      Knum <- ncol(pop_Q_vals)
      if (Knum == 2) { 
         k <- which.min( colSums(pop_Q_vals) ) 
         t_pop_Q_vals <- pop_Q_vals[,k]
         #t_pop_Q_vals <- as.matrix(t_pop_Q_vals, nrow=length(unique(pop)))
         t_pop_Q_vals <- cbind(names(t_pop_Q_vals), t_pop_Q_vals)
         colnames(t_pop_Q_vals) <- c("sites", paste("Qprops", 1:(Knum-1), sep="") )
      }

      if (Knum > 2) { 
         ikmax <- which.max( colSums(pop_Q_vals) ) 
         t_pop_Q_vals <- pop_Q_vals[,-ikmax]
         #t_pop_Q_vals <- as.matrix(t_pop_Q_vals, nrow=length(unique(pop)))
         t_pop_Q_vals <- cbind(rownames(t_pop_Q_vals), t_pop_Q_vals)
         colnames(t_pop_Q_vals) <- c("sites", paste("Qprops", 1:(Knum-1), sep="") )
      }
      return(t_pop_Q_vals)
   }


   trans_pop_Q <- transform_pop_Q(pop_Q_vals)

   pop_Q_file  <- paste(basedir,species,"/popgen/",treatment,"/lea/",species,"_",dataset,"_popQ.txt", sep="")
   write.table(trans_pop_Q, pop_Q_file, quote=FALSE, col.names=TRUE, row.names=FALSE)

   return(pop_Q_file)
}

