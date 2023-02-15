#' Import DArT (onerow) xls file and convert to csv 
#'
#' @param basedir  -- base directory for analyses 
#' @param species  -- a species name (name of directory within base directory) 
#' @param dataset  -- dataset name assigned by DArT (appears in csv file name) 
#' @param sheet    -- sheet number for singlerow dataset
#' @param fnum     -- number to append at end of filename
#' @return Name of a csv file, containing DArT genotype data that arrived in xls format. 
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' \dontrun{
#' dart_data            <- read_dart_onerow_csv(RandRbase,species,dataset)
#' } 

read_dart_xls2csv <- function(basedir,species,dataset, sheet=4, fnum=2) {

   require(readxl)                                                 
   datafile  <- paste(basedir,species,"/Report-",dataset,".xlsx",sep="")
   csvfile   <- paste(basedir,species,"/Report_",dataset,"_SNP_singlerow_", fnum,".csv",sep="")
   cat("\n")
   cat(" Reading data file:", datafile,"\n")
   
   ### Open and parse dart datafile
  
   x <- read_excel(datafile, sheet=sheet,  col_names=FALSE)
   y <- as.matrix(x)
   z <- gsub("," , "" , y)

   cat(" Writing data file:", csvfile,"\n")
   write.table(z,csvfile,quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")

   return(csvfile)
}


