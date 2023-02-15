#' Performs spatial krig of population ancestry coefficients
#' for prediction using R package gdm
#'
#' @param gdm_dir -- directory where gdm was performed
#' @param pixels_per_degree -- spatial resolution
#' @param buff    -- spatial buffer
#' @param krig_lambda       -- parameter for krig, if NULL, uses default
#' @return a list of rasters for use in predicting gdm on krigged ancestry coefficients
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' gdm_files <- dart2gdm(gms, etc)
#' }

krig_SNMF_gdm <- function(gdm_dir, pixels_per_degree=20, buff=0.5, krig_lambda=NULL, srast=NULL, Q=TRUE) {

   require(spam)
   require(fields); require(sp); require(raster)

   ppd    <- pixels_per_degree

   if (is.null(srast)) {
      srf     <- make_srast(gdm_dir, pixels_per_degree=pixels_per_degree, buff=buff, Q=Q)
   } else {
      srf     <- srast
   }

   sr <- raster(srf)

   gridlist <- list(seq( sr@extent@xmin, sr@extent@xmax, length=sr@ncols) , seq( sr@extent@ymin, sr@extent@ymax, length=sr@nrows))
   sgrid     <- make.surface.grid(gridlist)

   ### pairwise Fst
   results <- read.table( paste(gdm_dir,"/environ_Q_data.txt", sep=""), header=TRUE, sep=" " )
   Qind    <- which(substr(colnames(results),1,6)=="Qprops")

   if (length(Qind) >= 1) {
      cat("   krigging ", length(Qind), " Q variables ... \n" )
   } else {
      cat("   could not find any Q proportions to krig... stopping"); stop()
   }

   Qrasters <- list()

   for (qq in Qind) {

      qname <- colnames(results)[qq]
      X.u <- results[, 2:3]
      Friction <- results[, qq]

      fit <- Krig(X.u, Friction, lambda=krig_lambda)

      fit_over_grid   <- predict(fit, sgrid)
      fit_grid        <- as.surface(sgrid, fit_over_grid)

      mat <- t(fit_grid$z)
      mat <- mat[rev(1:nrow(mat)),]

      # Make a raster object for the localdiff surface
      rr <- raster(mat,
               xmn = min(fit_grid$x),
               xmx = max(fit_grid$x),
               ymn = min(fit_grid$y),
               ymx = max(fit_grid$y))

      rastfil  <-  paste(gdm_dir, "/", qname, "_model.tif",sep="")
      writeRaster(rr, rastfil, format="GTiff", overwrite=TRUE)

      Qrasters[[qname]] <- rastfil

   }
   return(Qrasters)
}
