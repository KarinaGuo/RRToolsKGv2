#' calls cdhit
#' 
#' @return cdhit output
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' ld_files <- dart2localdiff(gms, etc)
#' }

cluster_sequences_cdhit <- function(fa_fil, cl_fil, cutoff=cutoff, cdhdir) {

      fc <- file.create(cdh_fil <- tempfile("cdh"))

      cdhitbin  <- paste0(cdhdir,"/cd-hit")
      cdhitcall <- paste0(cdhitbin, " -i ", fa_fil, " -o ", cdh_fil, " -c ", cutoff)

      system(cdhitcall)
  
      cdh_out <- paste0(cdh_fil,".clstr")
      if ( file.exists(cdh_out) ) {

         cdh.df <- parse_cdhit_file(cdh_out)
         write.table(cdh.df,cl_fil,col.names=FALSE, row.names=FALSE,quote=FALSE,sep=",")
         unlist(cdh_out); unlist(cdh_fil); 

      } else {
         cat("cdhit outfile is missing... \n")
      }

      return(cl_fil)      
}



