#' Write a dart data list object to an R workspace file  
#'
#' run_sunder() receives dart data in list format, and reports a series of quality stats.
#' These are written to RandR/data/species/qual_stats/name*
#'
#' @param sm    -- list of SpaceMix data and params  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' run_SpaceMix(sm)
#
 

run_SpaceMix <- function(sm) {

   require(SpaceMix)

   sm_run <- run.spacemix.analysis(n.fast.reps=sm$nfr, fast.MCMC.ngen=sm$fng, fast.model.option=sm$fmo,
       long.model.option=sm$fmo, data.type=sm$dt, counts = sm$counts, sample.sizes = sm$sample_sizes,
       spatial.prior.X.coordinates=as.vector(sm$locations_x), spatial.prior.Y.coordinates=as.vector(sm$locations_y), round.earth=FALSE,
       long.run.initial.parameters = NULL, k=nrow(sm$counts), loci=ncol(sm$counts), ngen=200000, printfreq=100, samplefreq=100,
       mixing.diagn.freq = 100, savefreq=100, directory=sm$sm_dir, prefix=sm$prefix)
}



