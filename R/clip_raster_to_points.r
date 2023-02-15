#' Clip a raster to a convex hull surrounding a set of points
#'
#' @param basedir -- name of the base directory for R&R
#' @param species -- species name
#' @param dataset -- dataset name
#' @param treatment  
#' @param pop     -- a vector of population assignments
#' @param Ksel    -- the selected number of ancestral populations  
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' new_raster <- clip_raster_to_points(old_raster, location_file="locations.txt")
#' }

clip_raster_to_points <- function(raster_to_clip, location_file=NULL, gdm_dir=NULL, hull_res=150, convexity=0.85) {


   require(fields); require(sp); require(raster); require(INLA); require(rgeos)

   if (is.null(location_file)) {

      if (is.null(gdm_dir)) {
         cat("   No location file provided, and no gdm directory... stopping...\n"); stop()   
      }

      location_file   <- paste(gdm_dir, "/environ_data.txt",sep="")
       
   }

   location_info   <- read.table(location_file, header=TRUE, sep=" ")

   ld_hull         <- inla.nonconvex.hull(as.matrix(cbind(location_info$long, location_info$lat)), resolution = hull_res, convex = convexity)
   ld_poly         <- SpatialPolygons(list(Polygons(list(Polygon(ld_hull$loc)), ID=1)))
   ld_poly_buff    <- gBuffer(ld_poly,width=.1)

   clipped_raster <- mask(raster_to_clip, ld_poly_buff) 

   return(clipped_raster)
}
