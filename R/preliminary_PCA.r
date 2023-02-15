#' Write a dart data list object to an R workspace file  
#'
#' @param gl_file   -- location of an R file containing genlight object  [required]
#' @param basedir   -- Base dirctory  [required]
#' @param species   -- species name  [required]
#' @param dataset   -- dataset identifier  [required]
#' @param treatment -- treatment  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' preliminary_PCA(gl_file, basedir, species, dataset, treatment)
#
 

preliminary_PCA <- function(gl_file, data_list, basedir, species, dataset, treatment) {

   require(adegenet)

   if (file.exists(gl_file)) {

      cat("   Opening ", gl_file, "for PCA \n")
      load(gl_file)   

   } else {
      cat("   Error: The file ", gl_file, " does not exist. \n"); stop();
   }

   gl_dir    <- paste(basedir,species,"/popgen/",treatment,"/genlight", sep="")
   
   if(!dir.exists(gl_dir)) {
      cat("  The nominated directory ", gl_dir, " does not exist \n")
   } 

   pca_text   <- paste(gl_dir,"/",species,"_",dataset,"_PCA.txt",sep="")
   pca_file   <- paste(gl_dir,"/",species,"_",dataset,"_PCA.rda",sep="")
   pca_plot   <- paste(gl_dir,"/",species,"_",dataset,"_PCA.pdf",sep="")

   gl_pca <- glPca(dart_gl, nf=5, parallel=FALSE)

   pdf(file=pca_plot)
      scatter(gl_pca)
   dev.off()

   pca_table <- cbind(data_list$meta$site,data_list$meta$lat,data_list$meta$long,gl_pca$scores)
   colnames(pca_table)[1:3] <- c("site", "lat", "long")

   write.table(pca_table, file=pca_text, col.names=TRUE, row.names=TRUE, quote=FALSE, sep=",")
   
   save(gl_pca, file=pca_file)
}



