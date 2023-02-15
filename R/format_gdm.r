#' Prepares SNMF Q matrices to act as explanatory variables in GDM
#'
#' @param gdm_dir -- name of directory containing GDM data
#' @param Q       -- if TRUE, use Q matrix in GDM 
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_object <- format_gdm(gdm_dir, Q=FALSE)
#' }

format_gdm <- function(gdm_dir, Q=FALSE) {


   require(gdm)

   gdissim_file    <- paste(gdm_dir, "/allele_fst.txt",sep="")
   environ_file    <- paste(gdm_dir, "/environ_data.txt", sep="")

   if (Q == TRUE)  { environ_file    <- paste(gdm_dir, "/environ_Q_data.txt", sep="") }

   gdissim         <- read.table(gdissim_file, header=FALSE, sep=" ") 
   environ         <- read.table(environ_file, header=TRUE, sep=" ")

   id_gdissim      <- cbind(environ$sites, gdissim)
   colnames(id_gdissim)[1] <- "sites"

   gdm_spt <- formatsitepair(bioData=id_gdissim, bioFormat=3, predData=environ, siteColumn="sites", XColumn="long", YColumn="lat")

   gdm_in <- list(gdissim=gdissim, environ=environ, gdm_spt=gdm_spt)

   return(gdm_in)

}
