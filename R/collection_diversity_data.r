#' Calculates diversity (proportion of sites polymorphism) for
#' collections of sites
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms -- a dart data object [Required]
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param pop     -- a vector containing population assignments for each sample  
#' @return a list containing names of infiles for range expansion analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' dart_gi <- dart2rangexp(dms, meta.csv)
#' }

collection_diversity_data <- function(dms, basedir, species, dataset, pop, sites_single=TRUE, sites_paired=TRUE) {

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
   ind_NA_loci <- which( colSums(is.na(population_allele_stats$count)) > 0 )
   if ( length(ind_NA_loci) > 0 ) {
      cat("found ",  length(ind_NA_loci), "loci with no data for a population. Removing these loci \n")
      population_allele_stats$minor  <- population_allele_stats$minor[,-ind_NA_loci]
      population_allele_stats$count  <- population_allele_stats$count[,-ind_NA_loci]
      population_allele_stats$sample <- population_allele_stats$sample[,-ind_NA_loci]
      population_allele_stats$freq   <- population_allele_stats$freq[,-ind_NA_loci]
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

   cd_dir    <- paste(basedir,species,"/popgen/",treatment,"/CollecDiv", sep="")
   
   if(!dir.exists(cd_dir)) {
      cat("  collecdiv directory: ", cd_dir, " does not exist and is being created. \n")
      dir.create(cd_dir)
   } else {
      cat("  collecdiv directory: ", cd_dir, " already exists, content will be overwritten. \n")
   }

   counts       <- population_allele_stats$minor
   sample_sizes <- population_allele_stats$sample

   sm <- list(counts=counts, sample_sizes=sample_sizes)
   cd <- list(sm=sm, cd_dir=cd_dir, names=rownames(counts))

   AlleleCounts2CollecDiv <- function(sm, schemes) {

      ns <- nrow(schemes)
      np <- ncol(schemes)

      allele_counts <- mat.or.vec(ns,1)

      # loop through schemes and calculate allele counts
      for (i in 1:ns) {

         ps <- 0 
         for (j in 1:np) {
            if (schemes[i,j]==1) {
               if (ps == 0) {
                  scheme_counts     <- sm$counts[j,]
                  scheme_sampsizes  <- sm$sample_sizes[j,]
                  ps <- ps + 1
               }
               if (ps > 0) {
                  scheme_counts     <- scheme_counts + sm$counts[j,]
                  scheme_sampsizes  <- scheme_sampsizes + sm$sample_sizes[j,]
                  ps <- ps + 1
               }
            }
         }

         schemePolymorphic <- (scheme_counts!=0 & scheme_counts!=scheme_sampsizes) 
         schemePolymorphic[schemePolymorphic] <- 1
         schemePolyCount   <- sum(schemePolymorphic)              
         allele_counts[i]  <- schemePolyCount

      } # close loop through schemes

      return(allele_counts)
   }


   make_paired_collection_schemes <- function(sm) {
      n <- nrow(sm$counts)
      scs_pairs <- mat.or.vec(n*(n-1)/2,n)
      s <- 1
      for (i in 1:n) {
         for (j in 1:n) {
            if (i < j) {
               scs_pairs[s,c(i,j)] <- 1
               s <- s+1
            }
         } 
      }
      return(scs_pairs)
   }


   make_site_collection_schemes <- function(sm) {
      n <- nrow(sm$counts)
      scs_singles <- diag(n)
      return(scs_singles)
   }

   if (sites_single) {

      single_schemes   <- make_site_collection_schemes(sm)
      single_site_div  <- AlleleCounts2CollecDiv(sm, single_schemes)
      cd               <- c(cd, list(single_site_div=single_site_div))
   }

   if (sites_paired) {

      paired_schemes   <- make_paired_collection_schemes(sm)
      paired_site_div  <- AlleleCounts2CollecDiv(sm, paired_schemes)
      cd               <- c(cd, list(paired_site_div=paired_site_div))
   }

   return(cd)
}



