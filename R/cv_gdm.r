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

cv_gdm <- function(gdm_dir, geo=TRUE, Q=FALSE, krig_lambda=NULL, E=FALSE, erast=NULL, random_pairs=FALSE, number_random=100) {

   master_fst_file   <- paste(gdm_dir, "/allele_fst.txt",sep="")
   

   if (Q) { 
      master_env_file   <- paste(gdm_dir, "/environ_Q_data.txt",sep="") 
   } else {
      master_env_file   <- paste(gdm_dir, "/environ_data.txt",sep="")
   }

   # check files exist
   if(file.exists(master_fst_file) & file.exists(master_env_file)) {
      cat("")
   } else {
        cat("  One or both files needed for cross-validation of gdm are missing... \n")
   }

   # make a new directory
   CV_dir <- paste(gdm_dir,"/cv_pairs", sep="")
   if(!dir.exists(CV_dir)) {
      cat("  CV directory: ", CV_dir, " does not exist and is being created. \n")
      dir.create(CV_dir)
   } else {
      cat("  CV directory: ", CV_dir, " already exists, content will be overwritten. \n")
   }

   # read in data
   master_fst <- read.table(master_fst_file, header=FALSE, sep=" ")
   master_env <- read.table(master_env_file, header=TRUE, sep=" ")

   npop <- nrow(master_env)

   population_pairs <- mat.or.vec(npop*(npop-1)/2,2)
   s <- 1
   for (i in 1:npop) {
      for (j in 1:npop) {
         if (i < j) {
            population_pairs[s,1] <- i
            population_pairs[s,2] <- j
            s <- s+1
         }
      } 
   }

   if ( !random_pairs | number_random >= nrow(population_pairs) ) {
      cv_pairs <- population_pairs

   } else {

      rind <- sample(1:nrow(population_pairs))[1:number_random]
      cv_pairs <- population_pairs[ rind, ]
   }

   cv_results <- mat.or.vec(nrow(cv_pairs),2)
   colnames(cv_results) <- c("observed", "model")

   loc_vec <- NULL

   for (p in 1:nrow(cv_pairs)) {

      i_s1 <- cv_pairs[ p, 1 ]
      i_s2 <- cv_pairs[ p, 2 ]

      locations <- paste(master_env[i_s1,1],master_env[i_s2,1],sep="-")
      loc_vec   <- c(loc_vec, locations)
      #cv_results[p, 3] <- locations
      cv_results[p, 1] <- master_fst[i_s1, i_s2]
      
      p_dir <- paste(CV_dir,"/CV-",locations, sep="")
      dir.create(p_dir)

      i_pair <- c(i_s1, i_s2)
      p_ex_fst <- master_fst[-i_pair, -i_pair]  
      p_ex_env <- master_env[-i_pair, ]
  
      p_ex_fst_file <- paste(p_dir, "/allele_fst.txt", sep="")

      if (!Q) {
         p_ex_env_file <- paste(p_dir, "/environ_data.txt", sep="")
         qrast <- NULL 
      } else {
         p_ex_env_file <- paste(p_dir, "/environ_Q_data.txt", sep="")
      }

      write.table(p_ex_fst, p_ex_fst_file, col.names=FALSE, row.names=FALSE, quote=FALSE)
      write.table(p_ex_env, p_ex_env_file, col.names=TRUE, row.names=FALSE, quote=FALSE)

      if (Q) {
         # switched to krigcomps 26 Sept 2018
         # qrast <- krig_SNMF_gdm(p_dir)
         qrast <- krigcomps_SNMF_gdm(p_dir)
      }

      gdm_cv_p_ex      <- format_gdm(p_dir, Q=Q)
      gdm_format_file  <- paste(p_dir, "/gdm_ex_formatted.txt", sep="")      
      write.table(gdm_cv_p_ex$gdm_spt, gdm_format_file, col.names=TRUE, row.names=FALSE, quote=FALSE)

      mod_cv_p_ex      <- gdm(gdm_cv_p_ex$gdm_spt,geo=geo)

      #cat(p_dir, "\n")

      if (Q == TRUE & is.null(qrast) ) {
         cv_results[p, 2] <- NA
      } else {
         p_gdm            <- predict_gdm(p_dir, mod_cv_p_ex, srast=NULL, Q=Q, qrast=qrast, E=E, erast=erast, s1_point=as.numeric(master_env[i_s1,2:3]), s2_point=as.numeric(master_env[i_s2,2:3]), point_to_point=TRUE)
         cv_results[p, 2] <- p_gdm
      }

      # debug
      if(is.na(p_gdm)) {
         cat("NA value for CV prediction: ", locations, "  \n")
         s1_point=as.numeric(master_env[i_s1,2:3])
         s2_point=as.numeric(master_env[i_s2,2:3])
         save(p_dir, mod_cv_p_ex, Q, qrast, E, erast, s1_point, s2_point, file=paste(p_dir,"/prediction_debug.Rd",sep=""))
      }
      #unlink(p_dir, recursive=TRUE)
   }

   rownames(cv_results) <- loc_vec
   return(cv_results)
}

