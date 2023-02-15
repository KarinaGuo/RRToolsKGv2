#' Gets data for a list of climate data 
#' for a set of locations
#'
#' @param varlist -- A file with two columns, 'var' and 'file'
#'                   var lists the variable names
#'                   file is the name of the gtif 
#' @param climdir -- name of the directory containing gtif files 
#' @param longlat -- a set of long and lat info for data points
#'
#' @return -- climate variable data for each location
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' VarClimLoc <- get_climate_data(varlist=RRStandardClimLayers, climdir="/home/jason/RandR/clim", longlat=long_lat.txt)
#' }


get_climate_data <- function(varfile, climdir, longlatpoints) {
   
   require(raster)

   # checks
   if(!dir.exists(climdir)) {
      cat("The nominated climate directory does not exist \n"); stop();
   }

   if(!file.exists(varfile)) {
      cat("The file listing climate variables does not exist \n"); stop();
   }

   long_lat       <- longlatpoints

   clim_tifs      <- read.table(varfile, header=TRUE)$file
   colnames(long_lat) <- c("long", "lat") 
   clim_tifs_path <- paste(climdir, clim_tifs,sep="/")

   # open the climate layers in a stack, using raster
   s <- stack(clim_tifs_path)
   
   clim <- extract(s,long_lat)

   return(clim)
}


