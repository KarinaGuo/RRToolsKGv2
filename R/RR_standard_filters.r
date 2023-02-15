#' Read genotype filtering parameters from a file
#'
#' @param paramdir    -- name of directory containing a file with a standard parameter set [required]
#' @param param_file  -- name of the file assigning values to parameters [required]
#' @return            -- A list of parameters, which can be used in various data filtering functions
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' \dontrun{
#' p <- RR_standard_filters(RandRparam, param_file)
#' } 

RR_standard_filters <- function(paramdir, param_file) {

   file <- paste(paramdir, param_file, sep="")

   cat("\n")
   cat(" Reading parameters from file:", file,"\n")

   source(file)

   p <- list(threshold_missing_loci=p_threshold_missing_loci,
          min_repro=p_min_repro,
          max_missing_samples=p_max_missing_samples,
          seed=p_seed,               
          popassign=p_popassign)

   return(p)
}
