#' identifies clusters of loci base on a reference allele sequence 
#'
#' @param basedir   -- name of base directory [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param dataset   -- arbitrary name for the dataset [required]
#' @param cutoff    -- cutoff for clustering algorithm [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' write.dart.data(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 
cluster_allele_seqs_reference <- function(basedir, species, dataset, cutoff, cdbin) {

   

   base_fa <- basedir
   base_cl <- basedir

   fa_dir <- paste(base_fa, species, "/fasta",sep="")

   if(!dir.exists(fa_dir)) {
      cat("   sequeunce directory: ", fa_dir, " does not exist and is being created. \n")
      dir.create(fa_dir)
   } else {
      cat("   sequence directory: ", fa_dir, " already exists... content might be overwritten. \n")
   }

   cl_dir <- paste(base_cl, species, "/clust",sep="")

   if(!dir.exists(cl_dir)) {
      cat("   cluster directory: ", cl_dir, " does not exist and is being created. \n")
      dir.create(cl_dir)
   } else {
      cat("   cluster directory: ", cl_dir, " already exists... content might be overwritten. \n")
   }


   allele_tnam   <- paste(fa_dir,"/",species,"_",dataset,"_trimmed.tnam.clip.fasta",sep="")
   fasta_fil     <- allele_tnam

   cluster_tnam <- paste(cl_dir,"/",species,"_",dataset,"_allele.tnam.reference.clstrs",sep="")

   cat("   Commencing clustering of reference alleles with decipher \n")
   #cluster_sequences_decipher(fa_fil=fasta_fil,cl_fil=cluster_tnam, cutoff=0.10, method="complete",np=8)
   cluster_sequences_cdhit(fa_fil=fasta_fil,cl_fil=cluster_tnam, cutoff=cutoff, cdbin)

   cat("   Clusting complete... Matching to locus list file... \n")
   locus_list <- paste(fa_dir,"/",species,"_",dataset,"_locus_list.csv",sep="")

   lt <- read.table(locus_list, header=FALSE, sep=",")
   tc <- read.table(cluster_tnam, header=FALSE, sep=",")

   i_lt_tc <- match(lt[,2], tc[,1])

   combined <- data.frame(locus=lt[,1],cluster=tc[i_lt_tc,2])
   cluster_fil  <- paste(cl_dir,"/",species,"_",dataset,"_allele_ref_clust.csv",sep="")

   write.table(combined, cluster_fil, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=",")
   
   out <- list(fasta=fasta_fil, clusters=cluster_fil, cutoff=cutoff)

   date_line    <- date()
   species_line <- paste("Species: ", species, sep="")
   dataset_line    <- paste("Dataset: ", dataset, sep="")
   fasta_line    <- paste("Fasta file: ", fasta_fil, sep="")
   cluster_line  <- paste("Cluster file: ", cluster_fil, sep="")
   method_line  <- paste("Method: ", cdbin, sep="")
   param_line  <- paste("cutoff: ", cutoff, sep="")

   report_file <- file(paste(cl_dir,"/reference_cluster.log",sep=""))
   report <- c("cluster_allele_seqs_reference",date_line, species_line, dataset_line, fasta_line, cluster_line, method_line, param_line, "\n")
   writeLines(report, report_file)
   close(report_file)

   return(out)

}



