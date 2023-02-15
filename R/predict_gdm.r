#' Predict dissimilarity using a GDM 
#'
#' @param gdm_dir   -- the directory containing the gdm
#' @param gdm_model -- the gdm model to use in predictions
#' @param srast     -- a spatial raster for gdm prediction
#' @param Q         -- Q-matrix (ancestry coefficients) used in gdm  
#' @param qrast     -- if Q, a raster of Q matrix values
#' @param E         -- Environmental variables used in the model
#' @param erast     -- if E, a raster of environmental values
#' @param s1_point  -- a location, c(long, lat), for predicting dissimilarity
#' @param s2_point  -- a location, c(long, lat), for predicting dissimilarity
#' @param point_to_point -- if TRUE, model dissimilarity between s1 and s2 
#' @param point_to_all   -- if TRUE, model dissimilarity between s1 and all points in srast 
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' new_raster <- clip_raster_to_points(old_raster, location_file="locations.txt")
#' }


predict_gdm <- function(gdm_dir, gdm_model, srast=NULL, Q=FALSE, qrast=NULL, E=FALSE, erast=NULL, s1_point=NULL, s2_point=NULL, point_to_point=FALSE, point_to_all=FALSE, all_to_all=FALSE) {

   require(gdm)

   if (point_to_point) {
      if( is.null(s1_point) | is.null(s2_point) ) {
         cat("   For point to point prediction, s1 and s2 must be specified \n"); stop(); 
      } else {
         s1_ll   <- matrix(s1_point, nrow=1)
         s2_ll   <- matrix(s2_point, nrow=1)
         gdm_prd <- cbind(matrix(c(1), nrow=1), matrix(c(1), nrow=1), s1_ll, s2_ll)
         gdm_prd <- as.data.frame(gdm_prd)
         colnames(gdm_prd)  <- c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord")

         if(Q) {
            require(raster)
            qdata     <- stack(qrast)

            s1_Q_vals <- extract(qdata, matrix(s1_ll,ncol=2))
            s2_Q_vals <- extract(qdata, matrix(s2_ll,ncol=2))

            s1_Q      <- matrix( s1_Q_vals, ncol=length(s1_Q_vals))
            s2_Q      <- matrix( s2_Q_vals, ncol=length(s2_Q_vals))
            gdm_prd   <- cbind(gdm_prd, s1_Q, s2_Q)

            colnames(s1_Q) <- paste("s1.Q", 1:length(s1_Q_vals), sep="")
            colnames(s2_Q) <- paste("s2.Q", 1:length(s2_Q_vals), sep="")
         }

         if(E) {
            edata <- stack(erast)
            enames <- names(edata)
            s1_E  <- extract(edata, s1_ll)
            s2_E  <- extract(edata, s2_ll)

            colnames(s1_E) <- paste("s1.", enames, 1:ncol(s1_E), sep="")
            colnames(s2_E) <- paste("s2.", enames, 1:ncol(s2_E), sep="")

            gdm_prd   <- cbind(gdm_prd, s1_E, s2_E)
            
            # dbg
            # write.table(file=paste(gdm_dir,"/prd.csv",sep=""),gdm_prd,sep=",",quote=FALSE, row.names=TRUE,col.names=TRUE)
            #
         }

         # do prediction for point to point
         gdm_prediction  <- predict(gdm_model, gdm_prd)
         return(gdm_prediction)
      }

   }

   if (point_to_all) {
      if( is.null(s1_point) ) {
         cat("   For point to all prediction, s1 must be specified. \n"); stop(); 
      }

      if( !file.exists(srast) ) {
         cat("   For point to all prediction, drast must specified. \n"); stop(); 
      }

      sdata <- stack(srast)
      r.pts <- rasterToPoints(sdata, spatial=TRUE)

      r.pts@data <- data.frame(long=coordinates(r.pts)[,1],
                         lat=coordinates(r.pts)[,2], r.pts@data)                         

      s2_ll      <- r.pts@data[,1:2]
      null_dist  <- rep(1, nrow(s2_ll))
      null_wght  <- rep(1, nrow(s2_ll))

      s1_ll     <- cbind(rep(s1_point[1], nrow(s2_ll)), rep(s1_point[2], nrow(s2_ll)))

      gdm_prd   <- cbind(null_dist, null_wght, s1_ll, s2_ll)
      colnames(gdm_prd)  <- c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord")
      gdm_prd   <- as.data.frame(gdm_prd)

      if(Q) {
         qdata <- stack(qrast)
         s1_Q  <- extract(qdata, s1_ll)
         s2_Q  <- extract(qdata, s2_ll)

         colnames(s1_Q) <- paste("s1.Q", 1:ncol(s1_Q), sep="")
         colnames(s2_Q) <- paste("s2.Q", 1:ncol(s2_Q), sep="")

         gdm_prd   <- cbind(gdm_prd, s1_Q, s2_Q)

      }

     if(E) {
         edata <- stack(erast)
         enames <- names(edata)
         s1_E  <- extract(edata, s1_ll)
         s2_E  <- extract(edata, s2_ll)

         colnames(s1_E) <- paste("s1.", enames, 1:ncol(s1_E), sep="")
         colnames(s2_E) <- paste("s2.", enames, 1:ncol(s2_E), sep="")

         gdm_prd   <- cbind(gdm_prd, s1_E, s2_E)
 
         # in case environmental variable undefined in parts of srast
         gdm_prd <- gdm_prd[ rowSums(is.na(gdm_prd)) <= 0, ]
         s2_ll   <- gdm_prd[,5:6]
      }


      gdm_prediction  <- predict(gdm_model, gdm_prd)
      point2all_rast  <- raster(sdata)
      point2all_rast  <- rasterize(s2_ll,
                               point2all_rast,
                               field=gdm_prediction)

      return(point2all_rast)

   }

}

