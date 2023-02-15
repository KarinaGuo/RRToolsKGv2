#' Save information on DArT blast hits to eucalypt genome 
#'
#' save.eucalypt.genome.hits() saves information on DArT hits to eucalypt genome
#' This function is expected to be called during initial extraction of data
#' by the function read.dart.xls.onerow
#'
#' @param x         -- matrix of DArT data, supplied by DArT, as exctracted from excel file [required]
#' @param basedir  -- name of csv file containing the DArT data [required]
#' @param species  -- name of the species in 'GenuSpec' format [required]
#' @param dataset  -- dataset name assigned by DArT, appears in xlsx file name [required]
#' @return A statement about whether list consisting of a matrix of genotypes, locus names (CloneID), sample names, Reproducibility scores (RepAvg)
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com), based on template authored by Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' dart_data <- read.dart.onerow(datafile="Report-DFi16-2154_oneline.csv", topskip=6, nmetavar=38, nas="-")
#

 
save.eucalypt.genome.hits <- function(x, basedir, species, dataset) {

   eucalypt_aln_colnames <- c("CloneID","Chrom_Eucalyptus_v10_phytozome", "ChromPos_Eucalyptus_v10_phytozome", "AlnCnt_Eucalyptus_v10_phytozome", "AlnEvalue_Eucalyptus_v10_phytozome")
   
   column_indices <- which(colnames(x) %in% eucalypt_aln_colnames)

   if ( length(column_indices) != 5 ) {
      cat(" Warning: problem finding alignments to eucalypt genome.  \n")      
      cat(" Required headings not found in DArT file. Will not return hits. \n\n")      
      eucagenhits <- "none"
      return(eucagenhits)

   } else {

      eucagenhits <- x[, column_indices ]

      dir <- paste(basedir, species, "/genome",sep="")
      if(!dir.exists(dir)) {
         cat(" Directory: ", dir, " does not exist and is being created. \n")
         dir.create(dir)
      } else {
         cat(" Directory: ", dir, " already exists... content might be overwritten. \n")
      }

      file   <- paste(dir,"/",species,"_",dataset,"_aln_eucagen.rda",sep="")

      cat(" Hits eucalyptus genome written to: ", file, "\n\n")
      save(eucagenhits, file=file)
      return(file)

   }

}
