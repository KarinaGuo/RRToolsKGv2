#' Perform leave-n-out cross validation for gdm 
#'
#' @param gdm_dir   -- the directory containing the gdm
#' @param gdm_model -- the gdm model to use in predictions
#' @param n         -- number left out, n=1 or n=2
#' @param do_all    -- if n=2, do all pairs?
#' @param draws     -- if n=2, and do all is false, nominate a number of draws  
#' @param rekrig    -- if TRUE, rekrig Q martix 
#' @return file name 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' cv_gdm_dir <- cv_gdm(gdm_dir, gdm_model, )
#' }

loocv_gdm <- function(gdm_dir, krig_lambda=NULL, Q=FALSE) {

   master_fst_file   <- paste(gdm_dir, "/allele_fst.txt",sep="")
   master_env_file   <- paste(gdm_dir, "/environ_data.txt",sep="")

   # check files exist

   # make a new directory
   LOOCV_dir <- paste(gdm_dir,"/LOOCV", sep="")
   if(!dir.exists(LOOCV_dir)) {
      cat("  LOOCV directory: ", LOOCV_dir, " does not exist and is being created. \n")
      dir.create(LOOCV_dir)
   } else {
      cat("  LOOCV directory: ", LOOCV_dir, " already exists, content will be overwritten. \n")
   }

   # read in data

   master_fst <- read.table(master_fst_file, header=FALSE, sep=" ")
   master_env <- read.table(master_env_file, header=TRUE, sep=" ")

   if (Q)  { 
      master_env    <- paste(gdm_dir, "/environ_Q_data.txt", sep="") 
   }

   npop <- nrow(master_env)

   left_out_values <- mat.or.vec(npop,1)

   for (p in 1:npop) {
   
     long_lat_p <- matrix(master_env[p,2:3],nrow=1)
     # create a directory for the LOO analysis
     exdir <- paste(LOOCV_dir,"/p", p,sep="")
     if(!dir.exists(exdir)) {
        cat("  directory: ", exdir, " does not exist and is being created. \n")
        dir.create(exdir)
      } else {
         cat("  directory: ", exdir, " already exists, content will be overwritten. \n")
      }

      exfst <- master_fst[-p, -p]  
      exenv <- master_env[-p, ]
  
      exfst_file <- paste(exdir, "/allele_fst.txt", sep="")

      if (Q) {
         exenv_fil <- "/environ_Q_data.txt"
      } else {
         exenv_fil <- "/environ_data.txt"
      }

      exenv_file <- paste(exdir, env_fil, sep="")
      write.table(exfst, exfst_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(exenv, exenv_file, col.names=TRUE, row.names=FALSE, quote=FALSE)

      gdm_cv_pex     <-  format_gdm(exdir, Q=FALSE)
      mod_cv_pex     <-  gdm(gdm_cv_pex$gdm_spt,geo=TRUE)

      srast <- make_srast(gdm_dir=exdir)

      ex_out           <- mat.or.vec(npop,2)
      colnames(ex_out) <- c("observed","model")
      ex_out[,1]       <- master_fst[, p]
 
      for (pp in 1:npop) {
         if (pp != p) {
            long_lat_pp <- matrix(master_env[pp,2:3],nrow=1)
            p_gdm <- predict_gdm(exdir, mod_cv_pex, srast=NULL, s1_point=as.numeric(master_env[p,2:3]), s2_point=as.numeric(master_env[pp,2:3]), point_to_point=TRUE)
                     #predict_gdm(exdir, mod_cv_pex, srast=srast, s1_point=c(150,-35), s2_point=c(149,-28), point_to_point=TRUE)
            ex_out[pp,2] <- p_gdm 
            cat(p, " ", pp, " ", p_gdm, " \n")
         }
      }

#      call_localdiff(exdir)
#      krig_localdiff(exdir, pixels_per_degree=20, ld_dist=0.05, buffer=1.0, krig_lambda=krig_lambda)
#      rast_localdiff(exdir)

#      p_raster_file <- paste(exdir, "/ld_model.tif",sep="")
#      r <- raster(p_raster_file)

#      p_ld_left_out_val <- extract(r, lat_long_p)   
#      left_out_values[p] <- p_ld_left_out_val    

   }

 #  ld_inferred_vales  <- read.table(paste(ld_dir,"/outfst.txt",sep=""))
 #  LOOCV_results      <- cbind(ld_inferred_vales, left_out_values)
 #  LOOCV_results_file <- paste(LOOCV_dir, "/LOOCV_results.csv", sep="")
 #  colnames(LOOCV_results) <- c("pop", "long", "lat", "LocalDiff", "LOLocalDiff")
 #  write.table(LOOCV_results, LOOCV_results_file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep=",")

 #  LOOCV_results_figure <- paste(LOOCV_dir, "/LOOCV_results.pdf", sep="")
 #  pdf(LOOCV_results_figure)
 #  plot(LOOCV_results[,4],LOOCV_results[,5], xlab="localdiff value", ylab="left out localdiff estimated")
 #  abline(a=0, b=1)
 #  dev.off()

}

