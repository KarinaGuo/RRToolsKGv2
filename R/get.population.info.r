#' Returns information about populations, extracted from meta data
#' 
#' Input is a meta data object
#' 
#' @param m -- a meta data object [Required]
#' @param p -- the population assignments [Required]
#' @return a list containing locations, population names, population number
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

get.population.info <- function(m, p, method="average") {

   pop_list <- unique(p)
   num_pops <- length(pop_list)

   c <- 1
   lon_lat <- NULL
   for (i in 1:num_pops) {

      pop     <- pop_list[i]
      ind_pop <- which(p == pop)
      
      plats <- m$lat[ind_pop]
      plons <- m$long[ind_pop]

      if (method=="centroid"){
         pcent <- centroid(cbind(plons, plats))
         plon  <- as.numeric(pcent[1,1])
         plat  <- as.numeric(pcent[1,2])
      }

      if (method=="average"){
         plon <- mean(plons)
         plat <- mean(plats) 
      }

      p_lon_lat <- c(plon, plat)
      if (c == 1) {lon_lat <- p_lon_lat}
      if (c > 1)  {lon_lat <- rbind(lon_lat, p_lon_lat)}
      c <- c + 1
   }
   rownames(lon_lat) <- pop_list
   colnames(lon_lat) <- c("long", "lat")  

   pop_info <- list(lon_lat=lon_lat, names=pop_list, number=num_pops)
   return(pop_info)   

}



