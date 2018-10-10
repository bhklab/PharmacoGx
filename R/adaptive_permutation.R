############################################################################################################
######## Choosing sampling parameters for adaptive permutation testing                          ############
########                                                                                        ############
######## Ronglin Che and Chad C. Brown                                                          ############
########                                                                                        ############
######## The adaptive permutation (reference below) is a computationally efficient permutation  ############
######## procedure that is ideal for testing a large number of associations with a multiple     ############
######## comparison correction (e.g. Bonferroni), for example in genome-wide association        ############
######## studies.  This method has been shown to be virtually equivilent to standard            ############
######## permutation testing for p-values achieving family-wide significance, thereby           ############
######## preserving type I error rates.  In addition, in the case of testing the significance   ############
######## of between group variation (e.g. testing for differences between genotypes in GWASs)   ############
######## adaptive permutation has no loss of power when compared to ANOVA, when the assumption  ############
######## of normally distributed error terms holds.                                             ############
########                                                                                        ############
######## The functions below describe how to choose the permutation permameters (maximal        ############
######## number of permutations and maximal number of successes) in order to achieve the        ############
######## desired level of precision for adaptive permutation.                                   ############
########                                                                                        ############
######## Besag, Julian, and Peter Clifford. "Sequential monte carlo p-values." Biometrika 78.2  ############
######## (1991): 301-304.                                                                       ############
############################################################################################################


# This function choose_b determines the maximal number of 
# permutations for either adaptive or standard permutation.  
# It is a function of alpha (the p-value you would like to 
# estimate) and c (a desired precision level of p-value 
# estimation), where the standard error of the estimated 
# p-value is c*alpha.

choose_b <- function(alpha, c) {
  error <- alpha * c
  B <- alpha*(1 - alpha) / (error^2)
  return(B)
}

# This function choose_r determines the number of test 
# statistics that should be sampled in adaptive permutation 
# sampling such that p-values at the desired level of 
# significance (alpha) will be sampled with a 68% CI 
# (about 1 standard error) contained within a specified 
# level of precision (c).

choose_r <- function(alpha, c) {
  error <- alpha * c
  R <- 0
  foundR <- FALSE
  while(!foundR) {
    R <- R + 1
    brange <- qnbinom(c(0.1586553, 0.8413447), R, alpha)
    pvalRange <- R / (R + brange)
    diff <- max(abs(pvalRange - alpha))
    if(diff < error) {
      foundR <- TRUE
    }
  }
  return(R)
}

# For example, if the p-value threshold for significance 
# is 5*10^(-5), and the desired precision is 0.1 (so that 
# the standard error is 5*10^(-6), then parameters are 
# specified by the call: choose_b(alpha=5e-5, c=0.1) 
# which returns 1999900, and choose_r(alpha=5e-5, c=0.1) 
# which returns 121.