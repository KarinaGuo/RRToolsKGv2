#' Write dart locus sequences to a fasta file  
#'
#' @param x         -- content of dart data file  [required]
#' @param basedir   -- name of base directory [required]
#' @param species   -- name of the species in GenuSpec format [required]
#' @param name      -- arbitrary name for the dataset [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' write.dart.data(dart_data, RandRbase, "FicuRubi", "DFi16-2154")
#
 
write_allele_seq_fasta <- function(x, basedir, species, dataset, clip=5) {

   fa_dir <- paste(basedir, species, "/fasta",sep="")

   if(!dir.exists(fa_dir)) {
      cat("   sequeunce directory: ", fa_dir, " does not exist and is being created. \n")
      dir.create(fa_dir)
   } else {
      cat("   sequence directory: ", fa_dir, " already exists... content might be overwritten. \n")
   }

   locus_labels <- as.character(x$AlleleID)
   locus_names  <- as.character(lapply(strsplit(as.matrix(locus_labels),split="[|]"), function(x) x[1]))
   locus_pos    <- as.integer(x$SnpPosition)
   locus_SNP    <- as.character(x$SNP)
   locus_nuc    <- as.character(lapply(strsplit(as.matrix(locus_SNP),split="[:]"), function(x) x[2]))

   loc_labels   <- paste(locus_names,locus_pos,gsub(">","",locus_nuc),sep="_")
   #loc_labels   <- paste0("L",1:length(x$AlleleID))

   allele_file   <- paste(fa_dir,"/",species,"_",dataset,"_allele.fasta",sep="")
   trimmed_file  <- paste(fa_dir,"/",species,"_",dataset,"_trimmed.fasta",sep="")
   allele_t_file   <- paste(fa_dir,"/",species,"_",dataset,"_allele.tnam.fasta",sep="")
   trimmed_t_file  <- paste(fa_dir,"/",species,"_",dataset,"_trimmed.tnam.fasta",sep="")
   allele_tc_file   <- paste(fa_dir,"/",species,"_",dataset,"_allele.tnam.clip.fasta",sep="")
   trimmed_tc_file  <- paste(fa_dir,"/",species,"_",dataset,"_trimmed.tnam.clip.fasta",sep="")

   loclist_file  <- paste(fa_dir,"/",species,"_",dataset,"_locus_list.csv",sep="")

   a <- x$AlleleSequence
   t <- x$TrimmedSequence 

   cnum = clip + 1

   aclip <- substring(a, cnum)
   tclip <- substring(t, cnum)

   afat <- paste(">",loc_labels,"\n",a, sep="")
   tfat <- paste(">",loc_labels,"\n",t, sep="")

   afatc <- paste(">",loc_labels,"\n",aclip, sep="")
   tfatc <- paste(">",loc_labels,"\n",tclip, sep="")

   afa <- paste(">",x$AlleleID,"\n",a, sep="")
   tfa <- paste(">",x$AlleleID,"\n",t, sep="")

   write.table(afa, allele_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(tfa, trimmed_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(afat, allele_t_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(tfat, trimmed_t_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(afatc, allele_tc_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(tfatc, trimmed_tc_file, row.names=FALSE, col.names=FALSE, quote=FALSE)


   write.table(cbind(x$AlleleID, loc_labels), loclist_file, row.names=FALSE, col.names=FALSE, quote=FALSE,sep=",")

   return(allele_file)

}



