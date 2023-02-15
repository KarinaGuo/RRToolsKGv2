#' Writes a dart data object to a genepop input file
#' 
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms -- a dart data object [Required]
#' @return an object 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gp <- dart2genepop(dart_data, meta=FALSE)
#' }

dart2genepop <- function(dms, basedir, species, dataset, pop, maf_val=0.05, pop_miss_na=TRUE, exclude_if_single=TRUE) {

   treatment <- dms$treatment 
   if (dms$encoding == "altcount") {

      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genepop. \n")

   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }


   # global allele stats
   allele_frequency <- calculate.AF.summary(dms)
   ind_small_maf <- which( (allele_frequency$P < maf_val) | allele_frequency$P > (1-maf_val) )
   if ( length(ind_small_maf ) > 0 ) {
      cat("found ",  length(ind_small_maf), "loci with MAF smaller than ", maf_val," and removing these loci \n")
      dms <- remove.snps.from.dart.data(dms, ind_small_maf, input_as_names=FALSE)
      #gp_gt <- gp_gt[,-ind_small_maf]
   }

   # do some filtering
   population_allele_stats  <- calculate.population.allele.stats(dms, pop)

   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
      dms <- remove.snps.from.dart.data(dms, ind_NA_loci, input_as_names=FALSE)
      #gp_gt <- gp_gt[,-ind_NA_loci]
   }

   gp_gt <- dms$gt
   gp_gt[ gp_gt == 0 ] <- "0101"
   gp_gt[ gp_gt == 1 ] <- "0102"
   gp_gt[ gp_gt == 2 ] <- "0202"

   gp_gt[ is.na(gp_gt) ] <- "0000"



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

   gp_dir    <- paste(basedir,species,"/popgen/",treatment,"/genepop", sep="")
   
   if(!dir.exists(gp_dir)) {
      cat("  genepop directory: ", gp_dir, " does not exist and is being created. \n")
      dir.create(gp_dir)
   } else {
      cat("  genepop directory: ", gp_dir, " already exists, content will be overwritten. \n")
   }

   gp_file   <- paste(gp_dir,"/",species,"_",dataset,".gen",sep="")
 
   nS <- nrow(gp_gt); vS <- 1:nS; mS <- cbind(vS, rownames(gp_gt))
   nL <- ncol(gp_gt); vL <- paste("L", 1:nL, sep=""); mL <- cbind(vL, colnames(gp_gt))

   # write header
   sink(gp_file)
   cat(c("genepop file: ", species, "with MAF lt ", maf_val,"\n"))
   sink()

   # write locus line
   write(c(vL), ncolumns=(nL), file=gp_file, sep="\n", append=TRUE)
   #sink(gp_file, append = TRUE); cat(c("\n")); sink()

   poplist <- unique(pop)
   numpop  <- length(poplist)

   for (p in 1:numpop)
   {
      ithpop <- poplist[p]
      is     <- which(pop == ithpop)
      if (length(is) > 1) {
      sink(gp_file, append = TRUE); cat(c("POP\n")); sink()
       
      write.table(cbind(paste( rownames(gp_gt)[is]," ,",sep=" ") , matrix(gp_gt[is,],nrow=length(is))), file=gp_file, sep=" ",quote=FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
      } else {
         cat("   population ", poplist[p], "has 1 sample, ", rownames(gp_gt)[is], " excluded from gp file \n"  )
      }
   }

   # write data, pop by pop
   gp_fields <- list(gp_file=gp_file,treatment=treatment)   

   return(gp_file)

}
