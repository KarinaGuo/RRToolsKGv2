#' Provides an estimate of genetic diversity based on  
#' muliplicities of polymorphism (MOP) in DArT loci
#' 
#' Input is a DArT style genotype matrix
#' 
#' @param gt      -- genotype matrix [Required]
#' @param sample  -- sample size, in chromosomes
#' @param method  -- can be "diags", "mean", "mle", or "abc"
#' @return a list containing diversity estimates
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' div_calcs <- diversity_estimates_MOP(gt, method="mean")

diversity_estimates_MOP <- function(gt, method="mean") {

   # process gt for MOP data 
   require(stringr)
   loc_names <- str_split_fixed(as.character(colnames(gt)),"\\|",2)[,1]
   num_loc   <- length(unique(loc_names))
   num_seg   <- ncol(gt)
   mop_table <- table(table(loc_names))

   # process sample info
   samples <- nrow(gt)
   num_chr <- samples * 2

   # refine mop table
   m_vals <- as.numeric(names(mop_table))
   c_vals <- as.numeric(mop_table)

   m_max  <- max(m_vals)

   mop_zt <- mat.or.vec(m_max,2)
   for (i in 1:m_max) {
      mop_zt[i,1] <- i-1
      if ( i %in% m_vals ) {
         mop_zt[i,2] <- c_vals[ which(m_vals %in% i) ]
      } else {
         mop_zt[i,2] <- 0
      }
   }

   df_obs <- data.frame(mop_zt)
   colnames(df_obs) <- c("M", "F")

   if( m_max <= 1 | df_obs[1,1] != 0 ) {
      cat("Fatal error: there seems to be a problem with the MOP table"); stop()
   } else {
      cat("\n   MOP table found, proceeding to calculation of diagnostics\n")
   }

   # How well does poisson fit MOP
   wm_segsites      <- num_seg/num_loc-1   
   segsites_01      <- df_obs[ df_obs[,1]==1, 2] / df_obs[df_obs[,1]==0,2]

   #exp_mop          <- table(rpois( sum(df_obs[,2]) , wm_segsites))
   #df_exp           <- data.frame(exp_mop)
   #colnames(df_exp) <- c("M", "F")

   #df_contingency <- merge(df_exp, df_obs, by="M", all=TRUE)
   #df_contingency[ is.na(df_contingency) ] <- 0

   #contingency_exp_obs <- cbind(df_contingency[,"F.x"], df_contingency[,"F.y"])
   #colnames(contingency_exp_obs) <- c("expected", "observed") 
   #chisq_test_exp_obs  <- chisq.test(contingency_exp_obs, simulate.p.value=TRUE)

   out <- list() 
   #out$contingency <- contingency_exp_obs
   #out$fit_to_poisson <- chisq_test_exp_obs

   if(method=="mean") {
      out$est_ss_per_loc <- wm_segsites
      out$equiv_loc_num  <- num_seg / wm_segsites
   }

   if(method=="ratio") {
      out$est_ss_per_loc <- segsites_01
      out$equiv_loc_num  <- num_seg / segsites_01
   }

   out$sample_size   <- num_chr
   out$num_tags      <- num_loc
   out$num_snps      <- num_seg

   calculate_harmonic_number <- function(s) {
             h <- 0
             for (i in 1:(s-1)) {
                h <- h + 1/i
             }
             return(h)
   }

   out$harmonic_number <- calculate_harmonic_number(samples)


   return(out)
}





