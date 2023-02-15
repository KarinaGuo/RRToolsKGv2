#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param bd_par    -- list of BEDASSLE params  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_BEDASSLE(bd_par)
#
 

run_BEDASSLE <- function(bd) {

   require(BEDASSLE)
   bd_run <- MCMC(counts=bd$counts, sample_sizes=bd$sample_sizes, D=bd$D, E=bd$E, k=bd$k, loci=bd$loci, delta=bd$delta, aD_stp=bd$aD_stp, aE_stp=bd$aE_stp, a2_stp=bd$a2_stp, thetas_stp=bd$thetas_stp, mu_stp=bd$mu_stp, ngen=bd$ngen, printfreq=bd$printfreq, savefreq=bd$savefreq, samplefreq=bd$samplefreq, directory=bd$directory, prefix=bd$prefix) 
   

}



