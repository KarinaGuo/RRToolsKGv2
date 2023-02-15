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

krigcomps_SNMF_gdm <- function(gdm_dir, pixels_per_degree=20, buff=0.5, krig_lambda=NULL, srast=NULL, Q=TRUE, par_sph=100) {

   require(spam); require(compositions)
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

   results <- read.table( paste(gdm_dir,"/environ_Q_data.txt", sep=""), header=TRUE, sep=" " )
   Qind    <- which(substr(colnames(results),1,6)=="Qprops")

   if (length(Qind) >= 1) {
      cat("   krigging ", length(Qind), " Q variables ... \n" )
   } else {
      cat("   could not find any Q proportions to krig... stopping"); stop()
   }


   coords <- cbind(results$long,results$lat)
   colnames(coords) <- c("X", "Y")

   q_obs <- results[,Qind]
   q_rem <- 1 - rowSums(as.matrix(q_obs,nrow=length(Qind)))

   comp   <- acomp(cbind(q_obs, q_rem))

   lrv    <- logratioVariogram(comp, coords, nbins = 6)

   cat("Warning: setting global value of gpsph for eval statement \n")
   gpsph <<- par_sph
   vgModel <- CompLinModCoReg(~nugget()+sph(gpsph), comp)
   fit    <- vgmFit2lrv(lrv, vgModel, mode="log", iterlim = 2000) 

   if (fit$nlm$code == 3 | fit$nlm$code == 4 | fit$nlm$code == 5) {

      return( NULL )

   } else {

   coordsNew <- expand.grid(x = seq(from = sr@extent@xmin, to = sr@extent@xmax, length = sr@ncols), 
   y = seq(from = sr@extent@ymin, to = sr@extent@ymax, length = sr@nrows)) 

   krig <- compOKriging(comp, coords, coordsNew, fit$vg, err = FALSE)

   Qrasters <- list()

   for (qi in 1:length(Qind)) {

      qq <- Qind[qi]
      qname <- colnames(results)[qq]

      mat <- t(matrix(krig$Z[,qi],nrow=sr@ncols))
      mat <- mat[rev(1:nrow(mat)),]

# Make a raster object for the localdiff surface
      rr <- raster(mat,
               xmn = min(krig$X[,1]),
               xmx = max(krig$X[,1]),
               ymn = min(krig$X[,2]),
               ymx = max(krig$X[,2]))

      rastfil  <-  paste(gdm_dir, "/", species, "_", qname, ".tif",sep="")
      writeRaster(rr, rastfil, format="GTiff", overwrite=TRUE)

      Qrasters[[qname]] <- rastfil

   }

   return(Qrasters)

   }
}
