#'Fits curves of the form E = E_inf + (1 - E_inf)/(1 + (c/EC50)^HS) to dose-response data points (c, E) given by the user
#'and returns a vector containing estimates for HS, E_inf, and EC50.
#'
#'By default, logLogisticRegression uses an L-BFGS algorithm to generate the fit. However, if
#'this fails to converge to solution, logLogisticRegression samples lattice points throughout the parameter space.
#'It then uses the lattice point with minimal least-squares residual as an initial guess for the optimal parameters,
#'passes this guess to drm, and re-attempts the optimization. If this still fails, logLogisticRegression uses the
#'PatternSearch algorithm to fit a log-logistic curve to the data.
#'
#'@param conc [vector] is a vector of drug concentrations.
#'@param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#'drug concentrations whose logarithms are in the corresponding entries of the log_conc, where viability 0
#'indicates that all cells died, and viability 1 indicates that the drug had no effect on the cells.   
#'@param density [vector] is a vector of length 3 whose components are the numbers of lattice points per unit
#'length along the HS-, E_inf-, and base-10 logarithm of the EC50-dimensions of the parameter space, respectively. 
#'@param step [vector] is a vector of length 3 whose entries are the initial step sizes in the HS, E_inf, and
#'base-10 logarithm of the EC50 dimensions, respectively, for the PatternSearch algorithm.
#'@param precision is a positive real number such that when the ratio of current step size to initial step
#'size falls below it, the PatternSearch algorithm terminates. A smaller value will cause LogisticPatternSearch
#'to take longer to complete optimization, but will produce a more accurate estimate for the fitted parameters.
#'@param lower_bounds [vector] is a vector of length 3 whose entries are the lower bounds on the HS, E_inf,
#'and base-10 logarithm of the EC50 parameters, respectively.
#'@param upper_bounds [vector] is a vector of length 3 whose entries are the upper bounds on the HS, E_inf,
#'and base-10 logarithm of the EC50 parameters, respectively.
#'@param scale is a positive real number specifying the shape parameter of the Cauchy distribution.
#'@param Cauchy_flag [logical], if true, uses MLE under an assumption of Cauchy-distributed errors
#'instead of sum-of-squared-residuals as the objective function for assessing goodness-of-fit of
#'dose-response curves to the data.
#'@param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data,
#'and that log10(EC50) should be returned instead of EC50.
#'@param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#'than a percentage, and that E_inf should be returned as a decimal rather than a percentage.
#'@param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before 
#'curve-fitting is performed.
#'@param verbose [logical], if true, causes warnings thrown by the function to be printed.
#'@return A vector containing estimates for HS, E_inf, and EC50
#'@export
#'@importFrom stats optim
logLogisticRegression <- function(conc,
                                  viability,
                                  density = c(2, 10, 2),
                                  step = .5 / density,
                                  precision = 0.05,
                                  lower_bounds = c(0, 0, -6),
                                  upper_bounds = c(4, 1, 6),
                                  scale = 0.07,
                                  Cauchy_flag = FALSE,
                                  conc_as_log = FALSE,
                                  viability_as_pct = TRUE,
                                  trunc = TRUE,
                                  verbose = FALSE) {
  conc <- as.numeric(conc[!is.na(conc)])
  viability <- as.numeric(viability[!is.na(viability)])
  ii <- which(conc == 0)
  if(length(ii) > 0) {
    conc <- conc[-ii]
    viability <- viability[-ii]
  }
  
  
  
  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  if (prod(is.finite(conc)) != 1) {
    print(conc)
    stop("Concentration vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(viability)) != 1) {
    print(viability)
    stop("Viability vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(step)) != 1) {
    print(step)
    stop("Step vector contains elements which are not positive real numbers.")
  }
  
  if (prod(is.finite(precision)) != 1) {
    print(precision)
    stop("Precision value is not a real number.")
  }
  
  if (prod(is.finite(lower_bounds)) != 1) {
    print(lower_bounds)
    stop("Lower bounds vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(upper_bounds)) != 1) {
    print(upper_bounds)
    stop("Upper bounds vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(density)) != 1) {
    print(density)
    stop("Density vector contains elements which are not real numbers.")
  }
  
  if (is.finite(scale) == FALSE) {
    print(scale)
    stop("Scale is not a real number.")
  }
  
  if (is.logical(Cauchy_flag) == FALSE) {
    print(Cauchy_flag)
    stop("Cauchy flag is not a logical.")
  }
  
  if (is.logical(conc_as_log) == FALSE) {
    print(conc_as_log)
    stop("'conc_as_log' is not a logical.")
  }
  
  if (is.logical(viability_as_pct) == FALSE) {
    print(viability_as_pct)
    stop("'viability_as_pct' is not a logical.")
  }
  
  if (is.logical(trunc) == FALSE) {
    print(trunc)
    stop("'trunc' is not a logical.")
  }
  
  if (min(upper_bounds - lower_bounds) < 0) {
    print(rbind(lower_bounds, upper_bounds))
    stop("Upper bounds on parameters do not exceed lower bounds.")
  }
  
  if (length(conc) != length(viability)) {
    print(conc)
    print(viability)
    stop("Log concentration vector is not of same length as viability vector.")
  }
  
  if (min(viability) < 0) {
    if (verbose == TRUE) {
      warning("Warning: Negative viability data.")
    }
  }
  
  if (max(viability) > 100 / (1 + 99 * viability_as_pct)) {
    if (verbose == TRUE) {
      warning("Warning: Viability data exceeds negative control.")
    }
  }
  
  if (min(density) <= 0) {
    print(density)
    stop("Lattice point density vector contains negative values.")
  }
  
  if (precision <= 0) {
    print(precision)
    stop("Negative precision value.")
  }
  
  if (min(step) <= 0) {
    print(step)
    stop("Step vector contains nonpositive numbers.")
  }
  
  if (scale <= 0) {
    print(scale)
    stop("Scale parameter is a nonpositive number.")
  }
  
  if (conc_as_log == FALSE && min(conc) < 0) {
    print(conc)
    print(conc_as_log)
    stop("Negative concentrations encountered. Concentration data may be inappropriate, or 'conc_as_log' flag may be set incorrectly.")
  }
  
  if (viability_as_pct == TRUE && max(viability) < 5) {
    print(viability)
    print(viability_as_pct)
    if (verbose == TRUE) {
      warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
    }
  }
  
  if (viability_as_pct == FALSE && max(viability) > 5) {
    print(viability)
    print(viability_as_pct)
    if (verbose == TRUE) {
      warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
    }
  }
  
  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  if (conc_as_log == FALSE) {
    log_conc <- log10(conc)
  } else {
    log_conc <- conc
  }
  
  if (viability_as_pct == TRUE) {
    viability <- viability / 100
  }
  
  if (trunc == TRUE) {
    viability[which(viability < 0)] <- 0
    viability[which(viability > 1)] <- 1
  }
  
  #GENERATE INITIAL GUESS BY OBJECTIVE FUNCTION EVALUATION AT LATTICE POINTS
  sieve_guess <- .meshEval(log_conc,
                           viability,
                           lower_bounds = lower_bounds,
                           upper_bounds = upper_bounds,
                           density = c(2, 10, 2),
                           scale = scale,
                           Cauchy_flag = Cauchy_flag)
  sieve_guess_residual <- .residual(log_conc,
                                    viability,
                                    pars = sieve_guess,
                                    scale = scale,
                                    Cauchy_flag = Cauchy_flag)
  
  #ATTEMPT TO REFINE GUESS WITH L-BFGS OPTIMIZATION
  guess <- tryCatch(optim(par = sieve_guess,
                          fn = function(x) {.residual(log_conc,
                                                      viability,
                                                      pars = x,
                                                      scale = scale,
                                                      Cauchy_flag = Cauchy_flag)
                          },
                          lower = lower_bounds,
                          upper = upper_bounds,
                          method = "L-BFGS-B")[[1]],
                    error = function(e) {
                      c(NA, NA, NA)
                    })
  guess_residual <- .residual(log_conc,
                              viability,
                              pars = guess,
                              scale = scale,
                              Cauchy_flag = Cauchy_flag)
  
  #CHECK SUCCESS OF L-BFGS OPTIMIZAITON AND RE-OPTIMIZE WITH A PATTERN SEARCH IF NECESSARY
  if (prod(is.na(guess)) == 1 || guess_residual >= sieve_guess_residual) {
    guess <- sieve_guess
    guess_residual <- sieve_guess_residual
    span <- 1
    
    while (span > precision) {
      neighbours <- rbind(guess, guess, guess, guess, guess, guess)
      neighbour_residuals <- matrix(NA, nrow=1, ncol=6)
      neighbours[1, 1] <- pmin(neighbours[1, 1] + span * step[1], upper_bounds[1])
      neighbours[2, 1] <- pmax(neighbours[2, 1] - span * step[1], lower_bounds[1])
      neighbours[3, 2] <- pmin(neighbours[3, 2] + span * step[2], upper_bounds[2])
      neighbours[4, 2] <- pmax(neighbours[4, 2] - span * step[2], lower_bounds[2])
      neighbours[5, 3] <- pmin(neighbours[5, 3] + span * step[3], upper_bounds[3])
      neighbours[6, 3] <- pmax(neighbours[6, 3] - span * step[3], lower_bounds[3])
      
      for (i in 1:nrow(neighbours)) {
        neighbour_residuals[i] <- .residual(log_conc,
                                            viability,
                                            pars = neighbours[i, ],
                                            scale = scale,
                                            Cauchy_flag = Cauchy_flag)
      }
      
      if (min(neighbour_residuals) < guess_residual) {
        guess <- neighbours[which.min(neighbour_residuals), ]
        guess_residual <- min(neighbour_residuals)
      } else {
        span <- span / 2
      }
    }
  }
  return(list("HS" = guess[1],
              "E_inf" = ifelse(viability_as_pct, 100 * guess[2], guess[2]),
              "EC50" = ifelse(conc_as_log, guess[3], 10 ^ guess[3])))
}
