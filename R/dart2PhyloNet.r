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

dart2PhyloNet <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop) {

   #require(ape)

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

   phynet_dir    <- paste(basedir,species,"/popgen/",treatment,"/phynet", sep="")
   
   if(!dir.exists(phynet_dir)) {
      cat("  RE directory: ", phynet_dir, " does not exist and is being created. \n")
      dir.create(phynet_dir)
   } else {
      cat("  RE directory: ", phynet_dir, " already exists, content will be overwritten. \n")
   }

   phynet_nexus_file   <- paste(phynet_dir,"/",species,"_",dataset,".nex",sep="")


   if(!add_pop) {

      

   } else {

     poptext <- "BEGIN PHYLONET;\nMCMC_BiMarkers -cl 500000 -bl 200000 -sf 500 -diploid -op -varytheta -pp 2.0 -ee 2.0 -mr 1 -pl 4 -esptheta -ptheta 0.3\n-sd 12345678\n"

     #svdq_gt_out        <- matrix(rep("", nrow(gt)*ncol(gt) ), nrow=nrow(gt))
     sample_order       <- NULL
     up <- unique(pop)
     taxlist <- "-taxa ("
     taxmap  <- "-tm <"

     c <- 1
     for (i in 1:length(up)) {

         p  <- up[i]
         is <- which(pop == p)
         slen <- length(is)
         taxmap <- paste0(taxmap, p, ":")         
         
         cs <- 1
         for (s in is) {

            sampname <- rownames(gt)[s]
            taxmap <- paste0(taxmap, sampname)
            if (cs < length(is)) { taxmap <- paste0(taxmap, ",") }
            if (cs == length(is) & c != length(rownames(gt))) {taxmap <- paste0(taxmap, ";") }

            taxlist <- paste0(taxlist, sampname)
            if (c < length(rownames(gt))) { taxlist <- paste0(taxlist, ",") }
            if (c == length(rownames(gt))) {taxlist <- paste0(taxlist, ")") }

            c <- c + 1
            cs <- cs + 1
         }

         if (i == length(up)) {taxmap <- paste0(taxmap,">;\n\nEND;") }
     }
   } 

   nexhead <- "#NEXUS\nBegin data;\n"
   write(nexhead,file=phynet_nexus_file,append=FALSE)
   write(paste0("Dimensions ntax=",nrow(gt)," nchar=",ncol(gt),";"),file=phynet_nexus_file,append=TRUE)
   write("Format datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n",file=phynet_nexus_file,append=TRUE)

   for (i in 1:nrow(gt)) {
      write(paste0(rownames(gt)[i]," ", paste0(as.character(gt[i,]),collapse="")),file=phynet_nexus_file,append=TRUE)
   }

   write("\n;End;\n",file=phynet_nexus_file,append=TRUE)

   write(poptext,file=phynet_nexus_file,append=TRUE)
   write(taxlist,file=phynet_nexus_file,append=TRUE)
   write(taxmap,file=phynet_nexus_file,append=TRUE)

   return(phynet_nexus_file)

}
