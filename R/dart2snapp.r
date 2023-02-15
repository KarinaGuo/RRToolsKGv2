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
#' dart_gi <- dart2snapp(dart_data, meta.csv)
#' }

dart2snapp <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop) {

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

   snapp_dir    <- paste(basedir,species,"/popgen/",treatment,"/snapp", sep="")
   
   if(!dir.exists(snapp_dir)) {
      cat("  RE directory: ", snapp_dir, " does not exist and is being created. \n")
      dir.create(snapp_dir)
   } else {
      cat("  RE directory: ", snapp_dir, " already exists, content will be overwritten. \n")
   }

   snapp_nexus_file   <- paste(snapp_dir,"/",species,"_",dataset,".nex",sep="")

   if(!add_pop) {
      sample_names <- rownames(gt)
   } else {
      cat(   "   adding taxon names to samples \n")
      sample_names <- paste(pop, rownames(gt), sep="_")
   } 

   rownames(gt) <- sample_names
   write.nexus.data(gt, snapp_nexus_file)

   cat("   nexus file written... changing headers for beautii... \n")
   snf      <- scan(snapp_nexus_file, what="character", sep="\n")
   idt      <- grep("DATATYPE", snf)
   snf[idt] <- gsub("DNA", "STANDARD", snf[idt]) 
   snf[idt] <- gsub("INTERLEAVE", "SYMBOLS=\"012\" INTERLEAVE", snf[idt]) 

   cat("   overwriting nexus file with adjusted FORMAT line... \n")
   write(snf, snapp_nexus_file)

   return(snapp_nexus_file)

}
