#' Import DArT genotypes and store in a list 
#'
#' read.dart.xls.onerow() opens the onerow data format supplied by DArT
#' it is assumed this will be the fourth sheet in an xlsx file 
#' It skips the first n=topskip lines, and assumes the next line contains  
#' sample labels followed immediately by genotype data.
#'
#' It reads the genotype calls in dart format (encoding=dart) and
#' can covert these to counts of the alternate allele (encoding=altcount)
#'
#' The data matrix, locus names, metadata, sample names, are written to a dart data
#' object. This is simply a list containing some specified fields
#'
#' @param basedir  -- name of csv file containing the DArT data [required]
#' @param species  -- name of the species in 'GenuSpec' format [required]
#' @param dataset  -- dataset name assigned by DArT, appears in xlsx file name [required]
#' @param topskip  -- number of rows to skip before the header row (containing the specimen identities) [required]
#' @param nmetvar  -- number of columns containing the locus metadata (e.g. CloneID, Reproducibility) [required]
#' @param nas      -- missing data character [default "-"]
#' @param altcount -- if TRUE, change genotype encoding to the count of 'alternate' allele  [default TRUE]
#' @param euchits  -- if TRUE, gets information on BLAST hits in Eucalylptus genome  [default TRUE]
#' @return A DArT data list, genotype matrix, locus names (CloneID), sample names, Reproducibility scores (RepAvg)
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com), based on template authored by Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#' dart_data            <- read.dart.xls.onerow(RandRbase,species,dataset,topskip, nmetavar, euchits=TRUE)
#' } 

read.dart.xls.onerow <- function(basedir,species,dataset,topskip, nmetavar, nas="-", altcount=TRUE, euchits=FALSE) {

   require(readxl)
   datafile  <- paste(basedir,species,"/dart_raw/Report-",dataset,".xlsx",sep="")

   cat("\n")
   cat(" Reading data file:", datafile,"\n")

   ### Open and parse dart datafile
   # sheet 4 corresponds to onerow format
   x <- read_excel(datafile, sheet=4,  skip=topskip, col_names=TRUE, na=nas)
   
   ### Error checks
   if(any(names(x) == "CloneID")) {
   cat("  includes key variable CloneID\n")
   } else {
      cat(" Fatal Error: Dataset does not include key variable CloneID!\n"); stop()
   }
   if(any(names(x) == "RepAvg")) {
     cat("  includes key variable RepAvg\n\n")
   } else {
     cat(" Warning: Dataset does not include variable RepAvg! Check for errors\n\n"); stop()
   }

   NACloneID <- is.na(x$CloneID)
   if (any(NACloneID) ) {
     num_NACloneID <- length( which( NACloneID ) )
     cat(" Found ", num_NACloneID," missing CloneID values. Removing from data. \n\n")
     x <- x[ -which( NACloneID ), ]   
   } else {
     cat(" No missing CloneID values. \n\n")
   }
   
   if (!euchits) {
      cat(" Ignoring information on DArT locus alignment to Eucalyptus genome \n\n")
   } else {
      cat(" Saving information on DArT locus alignment to Eucalyptus genome \n")
      aln_save_status <- save.eucalypt.genome.hits(x, basedir, species, dataset)
   }

   # read some important data columns
   locus_labels <- as.character(x$CloneID)
   locus_names  <- as.character(lapply(strsplit(as.matrix(locus_labels),split="[|]"), function(x) x[1]))
   locus_repro  <- as.numeric(x$RepAvg)
   locus_calls  <- as.numeric(x$CallRate)
   locus_pos    <- as.integer(x$SnpPosition)
   locus_SNP    <- as.character(x$SNP)
   locus_nuc    <- as.character(lapply(strsplit(as.matrix(locus_SNP),split="[:]"), function(x) x[2]))

   num_snps <- nrow(x)
   num_loci <- length(unique(locus_names))
   num_cols <- ncol(x)
   num_samp <- num_cols - nmetavar

   # print some stats
   cat(" Initial data scan -- \n")
   cat("   Samples: ",num_samp," \n")
   cat("   SNPs: ", num_snps," \n")
   cat("   Loci: ", num_loci," \n\n")

   # place samples in rows, loci in columns

   tgt <- x[,(nmetavar+1):num_cols]
   gt <- as.matrix(t(tgt))

   # organize names
   sample_names <- colnames(x)[(nmetavar+1):num_cols]
   colnames(gt) <- as.vector(locus_labels)


   # make list
   cat(" Creating a DArT data list containing:                       \n")
   cat("             Genotypes                    --  $gt            \n") 
   cat("             Sample Names                 --  $sample_names  \n")
   cat("             Locus Names                  --  $locus_names   \n")
   cat("             Locus Reproducibility Scores --  $locus_repro   \n")
   cat("             Position of SNP in locus     --  $locus_pos     \n")
   cat("             Data filtering treatments    --  $treatment     \n")
   cat("             Position of SNP in locus     --  $locus_pos     \n")
   cat("             Method of data encoding gt   --  $encoding      \n")
   cat("             Nucleotides in this SNP      --  $locus_nuc     \n\n")

   treatment <- "raw"
   encoding  <- "DArT"
   dart_data <- list(gt=gt, sample_names=sample_names, locus_names=locus_names, locus_repro=locus_repro, locus_pos=locus_pos, locus_nuc=locus_nuc, encoding=encoding, treatment=treatment)

   if (altcount) {
       dart_data <- encode.dart2altcount(dart_data)
   } else {
       cat(" Warning: genotypes encoded in dart onerow format, 1=hom alt, 2=het \n\n")
   }

   return(dart_data)
}


