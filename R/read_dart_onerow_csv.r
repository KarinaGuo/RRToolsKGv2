#' Import DArT (onerow) genotype data and stores them in a list 
#'
#' It reads the genotype calls in dart format (encoding=dart) and
#' can covert these to counts of the alternate allele (encoding=altcount)
#'
#' @param basedir  -- base directory for analyses [required]
#' @param species  -- a species name (name of directory within base directory) [required]
#' @param dataset  -- dataset name assigned by DArT (appears in csv file name) [required]
#' @param misschar -- missing data character [default "-"]
#' @param altcount -- if TRUE, change genotype encoding to the count of 'alternate' allele  [default TRUE]
#' @param seq2fa   -- if TRUE, write allele sequences to a fasta file
#' @return A DArT data list, genotype matrix, locus names (by ID), sample names, Reproducibility scores (RepAvg)
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' \dontrun{
#' dart_data            <- read_dart_onerow_csv(RandRbase,species,dataset)
#' } 

read_dart_onerow_csv <- function(basedir,species,dataset, misschar="-", altcount=TRUE, seq2fa=FALSE, fnum=2) {

                                                 
   datafile  <- paste(basedir,species,"/Report_",dataset,"_SNP_singlerow_", fnum,".csv",sep="")

   cat("\n")
   cat(" Reading data file:", datafile,"\n")

   ### Open and parse dart datafile
   #x <- read.table(datafile,sep=",",na="-", stringsAsFactors=FALSE)
   x <- read.delim(datafile,sep=",",na="-", stringsAsFactors=FALSE, header=FALSE)
     
   rows_before_gt <- max(which(x[,1]=="*"))
   cols_before_gt <- max(which(x[1,]=="*"))

   colnames(x) <- x[rows_before_gt+1,]
   x <- x[-(1:(rows_before_gt+1)),]

   ### Error checks
   if(any(colnames(x) == "AlleleID")) {
   
   } else {
      
      cat("   Dataset does not include variable AlleleID... Check names... \n");

      if(any(colnames(x) == "CloneID")) {
         cat("   CloneID column found... proceeding with CloneID as names of loci. \n");
         colnames(x)[colnames(x)=="CloneID"] <- "AlleleID"
      } else {
         cat("   Locus names cannot be assigned. Error. \n"); stop()
      }
   }
   if(any(colnames(x) == "RepAvg")) {
     
   } else {
     cat("   Dataset does not include variable RepAvg... Check names... \n"); stop()
   }

   NAID <- is.na(x$AlleleID)
   if (any(NAID) ) {
     num_NAID <- length( which( NAID ) )
     cat(" Found ", num_NAID," missing AlleleID values. Removing from data. \n\n")
     x <- x[ -which( NAID ), ]   
   } 
   
   if (!seq2fa) {
      
   } else {
      cat("   Writing sequences to a fasta file \n")
      seq_fname <- write_allele_seq_fasta(x, basedir, species, dataset)
   }

   # read some important data columns
   locus_labels <- as.character(x$AlleleID)
   locus_names  <- as.character(lapply(strsplit(as.matrix(locus_labels),split="[|]"), function(x) x[1]))
   locus_repro  <- as.numeric(x$RepAvg)
   locus_calls  <- as.numeric(x$CallRate)
   locus_pos    <- as.integer(x$SnpPosition)
   locus_SNP    <- as.character(x$SNP)
   locus_nuc    <- as.character(lapply(strsplit(as.matrix(locus_SNP),split="[:]"), function(x) x[2]))

   num_snps <- nrow(x)
   num_loci <- length(unique(locus_names))
   num_cols <- ncol(x)
   num_samp <- num_cols - cols_before_gt

   # print some stats
   cat(" Initial data scan -- \n")
   cat("   Samples: ",num_samp," \n")
   cat("   SNPs: ", num_snps," \n")
   cat("   Loci: ", num_loci," \n\n")

   # place samples in rows, loci in columns

   tgt <- x[,(cols_before_gt+1):num_cols]
   gt  <- as.matrix(t(tgt))
   class(gt) <- "numeric"

   # organize names
   sample_names <- colnames(x)[(cols_before_gt+1):num_cols]
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


