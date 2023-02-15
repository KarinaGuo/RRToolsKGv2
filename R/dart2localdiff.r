#' Prepares input files that are ready
#' for analysis using software localdiff
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop     -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for localdiff analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' ld_files <- dart2localdiff(gms, etc)
#' }

dart2localdiff <- function(dms, basedir, species, dataset, pop) {

   # Step 1, get the genotypes ready
   treatment <- dms$treatment 
   if (dms$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   # Population allele stats
   population_allele_stats  <- calculate.population.allele.stats(dms, pop)
   population_spatial_dist  <- population.pw.spatial.dist(dms, pop)
   population_allele_fst    <- population.pw.Fst(dms, pop, basedir,species,dataset,maf_val = 0.1, miss_val = 0.1)

   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
      population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
      population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
      population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
      population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
   }

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

   ld_dir    <- paste(basedir,species,"/popgen/",treatment,"/localdiff", sep="")
   
   if(!dir.exists(ld_dir)) {
      cat("  localdiff directory: ", ld_dir, " does not exist and is being created. \n")
      dir.create(ld_dir)
   } else {
      cat("  localdiff directory: ", ld_dir, " already exists, content will be overwritten. \n")
   }

   #ld_object_file   <- paste(ld_dir,"/",species,"_",dataset,".rda",sep="")

   allele_fst_matrix      <- 1-population_allele_fst$Fst
   allele_cor_matrix      <- cor(t(population_allele_stats$minor))
   long_lat               <- population_spatial_dist$pop_info$lon_lat
   labels                 <- population_spatial_dist$pop_info$names

   distances_file         <- paste(ld_dir,"/spatial_dist.txt",sep="")
   allele_fst_file        <- paste(ld_dir,"/allele_fst.txt",sep="")
   allele_cor_file        <- paste(ld_dir,"/allele_cor.txt",sep="")
   long_lat_file          <- paste(ld_dir,"/long_lat.txt",sep="")
   label_file             <- paste(ld_dir,"/labels.txt",sep="")

 
   write.table(population_spatial_dist$S, distances_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")
   write.table(allele_fst_matrix, allele_fst_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")
   write.table(allele_cor_matrix, allele_cor_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")
   write.table(long_lat, long_lat_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")
   write.table(matrix(labels, nrow=1), label_file, col.names=FALSE, row.names=FALSE, quote=FALSE,sep=" ")

   allele_fst_pdf         <- paste(ld_dir,"/",species,"_allele_fst.pdf",sep="")
   pdf(file=allele_fst_pdf)
   plot(  as.vector(population_spatial_dist$S) , as.vector(allele_fst_matrix), xlab="distance", ylab=paste(species, " Fst",sep="") )
   dev.off()

   classic_ibd_pdf         <- paste(ld_dir,"/",species,"_ibd_plot.pdf",sep="")
   pdf(file=classic_ibd_pdf)
   plot(  as.vector(population_spatial_dist$S) , as.vector(population_allele_fst$Fst) / as.vector(1-population_allele_fst$Fst), xlab="distance", ylab=paste(species, " Fst/(1-Fst)",sep="") )
   dev.off()

   allele_cor_pdf         <- paste(ld_dir,"/",species,"_allele_cor.pdf",sep="")
   pdf(file=allele_cor_pdf)
   plot(  as.vector(population_spatial_dist$S) , as.vector(allele_cor_matrix), xlab="distance", ylab=paste(species, " Cor",sep="") )
   dev.off()


   return(ld_dir)

}
