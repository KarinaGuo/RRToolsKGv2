#' Generates the JSFS for two nominated populations 
#'
#' Input is a dart data object with altcount genotype encoding
#' (0=hom ref, 1=het, 2=hom alt).
#' 
#' @param dms      -- a dart data object [Required]
#' @param basedir  -- name of the base directory for R&R
#' @param species  -- species name
#' @param dataset  -- dataset name
#' @param samples  -- list containing indices for the populations for JSFS  
#' @return a list containing names of infiles for localdiff analysis
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' ld_files <- dart2localdiff(gms, etc)
#' }

dart2JSFS <- function(dms, basedir, species, samples) {

   pboth <- c(samples$p1, samples$p2)   

   num_p1 <- length(samples$p1)
   num_p2 <- length(samples$p2)

   # process genetic data
   gt      <- dms$gt
   num_loc <- ncol(gt)
   num_sam <- nrow(gt)

   # get the genotypes ready
   treatment <- dms$treatment 
   if (dms$encoding == "altcount") {
      cat(" Dart data object for ", dataset, "in species", species, "\n")
      cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
   } else {
      cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
   }

   ind_fixed <- which(  colSums(gt[ pboth,  ]) == ((num_p1 + num_p2)*2) | colSums(gt[ pboth,  ]) == 0 )   
   ind_NA    <- which(  is.na(colSums(gt[ pboth,  ])) )
   ind_rm    <- union(ind_fixed, ind_NA)

   gt_SNP    <- gt[ pboth, -ind_rm]
   p1        <- which( (1:num_sam)[c(samples$p1,samples$p2)] %in% samples$p1)   
   p2        <- which( (1:num_sam)[c(samples$p1,samples$p2)] %in% samples$p2)

   SFS  <- mat.or.vec(1, (num_p1+num_p2)*2 + 1)
   p1SFS  <- mat.or.vec(1, num_p1*2 + 1)
   p2SFS  <- mat.or.vec(1, num_p2*2 + 1)
   JSFS   <- mat.or.vec(num_p1*2 + 1, num_p2*2 + 1)

   for ( i in 1:ncol(gt_SNP) ) {

       p1_count <- sum(gt_SNP[ p1, i]) 
       p2_count <- sum(gt_SNP[ p2, i]) 
       pb_count <- sum(gt_SNP[ c(p1, p2), i]) 

       if ( (p1_count + p2_count) > length(pboth)  ) {
           p1c <- 2*num_p1 - p1_count + 1
           p2c <- 2*num_p2 - p2_count + 1
           pbc <- 2*length(pboth) - pb_count + 1

       } else {
          p1c <- p1_count + 1
          p2c <- p2_count + 1
          pbc <- pb_count + 1
       }

       SFS[pbc] <- SFS[pbc] + 1
       p1SFS[p1c] <- p1SFS[p1c] + 1
       p2SFS[p2c] <- p2SFS[p2c] + 1
       JSFS[ p1c, p2c ] <- JSFS[ p1c, p2c ] + 1 
   }
   private <- c(sum(JSFS[1,]), sum(JSFS[,1]), sum(JSFS[-1,-1]))

   SFS <- list(SFS=SFS, p1SFS=p1SFS, p2SFS=p2SFS, JSFS=JSFS, indp1=samples$p1, indp2=samples$p2, private=private)

   return(SFS)

}
