#' Returns a dart data object in appropriate files for analysis 
#' with the Peter and Slatkin method of detecting a range expansion 
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dart_data -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2svdq(dart_data, meta.csv)
#' }

dart2svdq <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop) {

   require(ape)

   # Step 1, get the genotypes ready
   treatment <- dart_data$treatment 
   if (dart_data$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # genotypes -- change char for missing
   gt <- dart_data$gt
   gt[ is.na(gt) ] <- "?"

   svdq_2row_gt     <- matrix(rep("", nrow(gt)*2*ncol(gt) ), nrow=nrow(gt)*2)
   svdq_2row_names  <- matrix(rep("", nrow(gt)*2 ), nrow=nrow(gt)*2)

   for (i in 1:ncol(gt) ) {

      ref_allele <- substr(dart_data$locus_nuc[i],1,1)
      alt_allele <- substr(dart_data$locus_nuc[i],3,3)

      for ( j in 1:nrow(gt) ) {

           g <- gt[j,i]
           a1 <- (j-1)*2 + 1
           a2 <- a1 + 1

           if (g == 0)   { svdq_2row_gt[ a1,i ] <- ref_allele; svdq_2row_gt[ a2,i ] <- ref_allele  }
           if (g == 1)   { svdq_2row_gt[ a1,i ] <- ref_allele; svdq_2row_gt[ a2,i ] <- alt_allele  }
           if (g == 2)   { svdq_2row_gt[ a1,i ] <- alt_allele; svdq_2row_gt[ a2,i ] <- alt_allele  }
           if (g == "?") { svdq_2row_gt[ a1,i ] <- "?"; svdq_2row_gt[ a2,i ] <- "?"  }

      }

   }


   for ( k in 1:nrow(gt) ) {

        a1 <- (k-1)*2 + 1
        a2 <- a1 + 1

        svdq_2row_names[a1,1] <- paste0(rownames(gt)[k],"_1")
        svdq_2row_names[a2,1] <- paste0(rownames(gt)[k],"_2")

   }

   rownames(svdq_2row_gt) <- svdq_2row_names

   # make directory, write files 
   dir <- paste(basedir, species, "/popgen",sep="")
   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
   }

   dir <- paste(basedir, species, "/popgen/",treatment,sep="")

   if(!dir.exists(dir)) {
      cat("  Directory: ", dir, " does not exist and is being created. \n")
      dir.create(dir)
   } else {
      cat("  Directory: ", dir, " already exists...  \n")
   }

   svdq_dir    <- paste(basedir,species,"/popgen/",treatment,"/svdq", sep="")
   
   if(!dir.exists(svdq_dir)) {
      cat("  RE directory: ", svdq_dir, " does not exist and is being created. \n")
      dir.create(svdq_dir)
   } else {
      cat("  RE directory: ", svdq_dir, " already exists, content will be overwritten. \n")
   }

   svdq_nexus_file   <- paste(svdq_dir,"/",species,"_",dataset,".nex",sep="")

#   if(!add_pop) {
#      sample_names <- rownames(gt)
#   } else {
#      cat(   "   adding taxon names to samples \n")
#      sample_names <- paste(rownames(gt), pop, sep="_")
#   } 

   write.nexus.data(svdq_2row_gt, svdq_nexus_file)

   cat("   nexus file written... changing headers for beautii... \n")
   snf      <- scan(svdq_nexus_file, what="character", sep="\n")
   idt      <- grep("DATATYPE", snf)
   snf[idt] <- gsub("DNA", "STANDARD", snf[idt]) 
   snf[idt] <- gsub("INTERLEAVE", "SYMBOLS=\"012\" INTERLEAVE", snf[idt]) 

   cat("   overwriting nexus file with adjusted FORMAT line... \n")
   write(snf, svdq_nexus_file)

   return(svdq_nexus_file)

}
