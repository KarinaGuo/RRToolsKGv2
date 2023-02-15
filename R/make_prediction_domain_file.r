#' Clip a raster to a convex hull surrounding a set of points
#'
#' @param species       -- species name
#' @param location_file -- file containing longs and lats
#' @param gdm_dir       -- directory used to fit GDM  
#' @param hull_res      -- resolution of convex hull
#' @param convexity     -- convexity of hull around sampled points  
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' new_raster <- make_prediction_domain_file(old_raster, location_file="locations.txt")
#' }

make_prediction_domain_file <- function(species, out_dir, location_file=NULL, gdm_dir=NULL, hull_res=150, convexity=0.85, buff_width=0.1) {


   require(fields); require(sp); require(raster); require(INLA); require(rgeos); require(rgdal)

   if (is.null(location_file)) {

      if (is.null(gdm_dir)) {
         cat("   No location file provided, and no gdm directory... stopping...\n"); stop()   
      }

      location_file   <- paste(gdm_dir, "/environ_data.txt",sep="")
       
   }

   location_info   <- read.table(location_file, header=TRUE, sep=" ")

   ld_hull         <- inla.nonconvex.hull(as.matrix(cbind(location_info$long, location_info$lat)), resolution = hull_res, convex = convexity)
   ld_poly         <- SpatialPolygons(list(Polygons(list(Polygon(ld_hull$loc)), ID=1)))
   ld_poly_buff    <- gBuffer(ld_poly,width=buff_width)

   sp_df           <- data.frame( ID=1:length(ld_poly_buff))
   spdf            <- SpatialPolygonsDataFrame(ld_poly_buff, sp_df, match.ID=FALSE)
   vector_file     <- paste(species,"_domain",sep="")
   writeOGR(spdf,paste(out_dir,vector_file,sep=""), "GeoJSON", driver="GeoJSON", check_exists=TRUE, overwrite_layer=TRUE)
                                                      #"GeoJSON"
                                                      # output as WGS84
   return(ld_poly_buff)
}
