#' Make a plot of admixture coefficients from an sNMF analysis
#'
#' plot_SNMF() information about an sNMF run, opens the project
#' A pdf file of admixture coefficients is made for a nominated K, R
#'
#' @param basedir   -- Base dirctory   [required]
#' @param species   -- species name    [required]
#' @param dataset   -- dataset identifier  [required]
#' @param treatment -- treatment       [required]
#' @param K         -- sNMF K value    [required]
#' @param r         -- sNMF replicate  [required]
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)

plot_SNMF <- function(basedir, species, dataset, treatment, Kvals=2:6, Rvals=1:3) {

   require(LEA)

   f <- paste(basedir, species, "/popgen/",treatment,"/lea/",species,"_",dataset,".snmfProject", sep="")


   if (file.exists(f)) {
      cat("   plotting sNMF Q values for ", species, "\n") 
   } else {
      cat("   Error: the sNMF project", f," does not exist \n"); stop(); 
   }

   s <- load.snmfProject(f)

   c <- c("blue", "red", "yellow", "green", "light blue", "orange")


   for (K in Kvals) {
      for (R in Rvals) {

         p <- paste(basedir, species, "/popgen/",treatment,"/lea/",species,"_",dataset,"_K", K, "_R", R, ".pdf", sep="")
         q <- Q(s,K,R)
         cK <- c[1:K]

         pdf(file=p)
            barplot(t(q), col=cK, border = NA, space = 0, xlab = "Sample", ylab = "Q")
          dev.off()
      }
   }
   return(p)
}



