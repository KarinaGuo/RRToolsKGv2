calculate.I.p1.p2 <- function(p1,p2){

   q1    <- 1 - p1
   q2    <- 1 - p2

   pmean <- (1/2) * ( p1 + p2 )

   A1 <- p1 * log(  p1 / pmean   )
   a1 <- q1 * log(  q1 / ( 1 - pmean )   )
   A2 <- p2 * log(  p2 / pmean   )
   a2 <- q2 * log(  q2 / ( 1 - pmean )  )

   I <- (1/2) * (A1 + a1 + A2 + a2)
 
   return(I)

}
