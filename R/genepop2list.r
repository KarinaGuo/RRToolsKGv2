genepop2list <- function(filename, nloc){

   require(stringr)
   lines <- readLines(filename)
   ss <- substring(lines,1,4)
   samps <- lines[ which( ss != "This" & ss != "pop " & ss!= "locu") ]

   samps_mat <- str_split_fixed(samps," ",(nloc+3))
   x <- samps_mat[,1]
   y <- samps_mat[,2]
   p <- paste("p_",x,"_",y,sep="")
   gt_gp <- samps_mat[ , 4:ncol(samps_mat)]
   gt_gp <- substring(gt_gp,1,6)

   gt <- as.matrix(gt_gp)
   gt[ which(gt_gp == "001001")  ] <- 0
   gt[ which(gt_gp == "001002")  ] <- 1
   gt[ which(gt_gp == "002001")  ] <- 1
   gt[ which(gt_gp == "002002")  ] <- 2

   class(gt) <- "numeric"

   locus_names  <- paste("L", 1:nloc, sep="")
   sample_names <- paste("S", 1:nrow(gt), sep="")
   colnames(gt) <- locus_names
   rownames(gt) <- sample_names

   lat      <- x
   long     <- y
   class(lat) <- "numeric"
   class(long) <- "numeric"

   site     <- p
   analyses  <- p 
   meta <- list(lat=lat, long=long, site=site, analyses=analyses)

   data <- list(gt=gt, treatment="simulated", encoding="altcount", sample_names=sample_names, locus_names=locus_names, meta=meta)

   return(data)
}
