#' Creates a raster with chosen resolution and buffers, and arbitrary values
#' Used to create a raster to predict a gdm model onto
#'
#' @param gdm_dir -- directory where gdm was performed
#' @param pixels_per_degree -- spatial resolution
#' @param buff    -- spatial buffer
#' @return a list of rasters for use in predicting gdm on krigged ancestry coefficients
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- dart2gdm(gms, etc)
#' }

make_srast <- function(gdm_dir, pixels_per_degree=20, buff=0.5, return_rast=TRUE, Q=FALSE) {

   require(spam)
   require(fields); require(sp); require(raster)

   ppd    <- pixels_per_degree

   location_file   <- paste(gdm_dir, "/environ_data.txt",sep="")

   if (Q) {
      location_file   <- paste(gdm_dir, "/environ_Q_data.txt",sep="")
   }
   location_info   <- read.table(location_file, header=TRUE, sep=" ")
   lat_long        <- as.matrix(cbind(location_info$long, location_info$lat))
   
   minlon <- min(lat_long[, 1])
   maxlon <- max(lat_long[, 1])

   lon_deg_buff <- (maxlon-minlon)*buff
   minlonbuff   <- minlon - lon_deg_buff
   maxlonbuff   <- maxlon + lon_deg_buff

   roundminlonbuff <- floor(minlonbuff*ppd)/ppd
   roundmaxlonbuff <- ceiling(maxlonbuff*ppd)/ppd

   minlat <- min(lat_long[, 2])
   maxlat <- max(lat_long[, 2])

   lat_deg_buff <- (maxlat-minlat)*buff
   minlatbuff <- minlat - lat_deg_buff
   maxlatbuff <- maxlat + lat_deg_buff

   roundminlatbuff <- floor(minlatbuff*ppd)/ppd
   roundmaxlatbuff <- ceiling(maxlatbuff*ppd)/ppd

   np_lon <- (roundmaxlonbuff-roundminlonbuff)*ppd
   np_lat <- (roundmaxlatbuff-roundminlatbuff)*ppd

   if (return_rast) {
      r <- raster(ncol=np_lon, nrow=np_lat, xmx=roundmaxlonbuff, xmn=roundminlonbuff, ymn=roundminlatbuff, ymx=roundmaxlatbuff)
      values(r) <- 1
      srf  <-  paste(gdm_dir, "/srast.tif",sep="")
      writeRaster(r, srf, format="GTiff", overwrite=TRUE)
      return(srf)
   } else {
      gridinfo <- list(np_lon=np_lon, np_lat=np_lat, xmx=roundmaxlonbuff, xmn=roundminlonbuff, ymn=roundminlatbuff, ymx=roundmaxlatbuff)
      return(gridinfo)
   }
}
