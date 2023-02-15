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
#' dart_gi <- dart2genalex(dart_data, meta.csv)
#' }

dart2genalex <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop=NULL, numeric_gt=FALSE) {

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

   genalex_2col_gt     <- matrix(rep("", nrow(gt)*2*ncol(gt) ), ncol=ncol(gt)*2)
   genalex_2col_loc_names  <- matrix(rep("", ncol(gt)*2 ), ncol=ncol(gt)*2)

   for (i in 1:ncol(gt) ) {

      ref_allele <- substr(dart_data$locus_nuc[i],1,1)
      alt_allele <- substr(dart_data$locus_nuc[i],3,3)

      for ( j in 1:nrow(gt) ) {

           g <- gt[j,i]
           a1 <- (i-1)*2 + 1
           a2 <- a1 + 1

           if (g == 0)   { genalex_2col_gt[ j,a1 ] <- ref_allele; genalex_2col_gt[ j,a2 ] <- ref_allele  }
           if (g == 1)   { genalex_2col_gt[ j,a1 ] <- ref_allele; genalex_2col_gt[ j,a2 ] <- alt_allele  }
           if (g == 2)   { genalex_2col_gt[ j,a1 ] <- alt_allele; genalex_2col_gt[ j,a2 ] <- alt_allele  }
           if (g == "?") { genalex_2col_gt[ j,a1 ] <- "N"; genalex_2col_gt[ j,a2 ] <- "N"  }

           genalex_2col_loc_names[a1] <- colnames(gt)[i]

      }

   }

   colnames(genalex_2col_gt) <- genalex_2col_loc_names
   
   if(!add_pop) {

      rownames(genalex_2col_gt) <- rownames(gt)

   } else {

      rownames(genalex_2col_gt) <- paste(rownames(gt),pop,sep=",")

   }


   if(numeric_gt) {

      genalex_2col_gt[ genalex_2col_gt == "A" ] <- "1"
      genalex_2col_gt[ genalex_2col_gt == "C" ] <- "2"
      genalex_2col_gt[ genalex_2col_gt == "G" ] <- "3"
      genalex_2col_gt[ genalex_2col_gt == "T" ] <- "4"
      genalex_2col_gt[ genalex_2col_gt == "N" ] <- "0"

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

   gax_dir    <- paste(basedir,species,"/popgen/",treatment,"/genalex", sep="")
   
   if(!dir.exists(gax_dir)) {
      cat("  RE directory: ", gax_dir, " does not exist and is being created. \n")
      dir.create(gax_dir)
   } else {
      cat("  RE directory: ", gax_dir, " already exists, content will be overwritten. \n")
   }

   gax_file   <- paste(gax_dir,"/",species,"_",dataset,".csv",sep="")

   write.table(genalex_2col_gt, gax_file, quote=FALSE, col.names=TRUE, row.names=TRUE, sep=",")


   return(gax_file)

}
