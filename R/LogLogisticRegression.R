#' Fits curves of the form E = E_inf + (1 - E_inf)/(1 + (c/EC50)^HS) to dose-response data points (c, E) given by the user
#' and returns a vector containing estimates for HS, E_inf, and EC50.
#' 
#' By default, logLogisticRegression uses an L-BFGS algorithm to generate the fit. However, if
#' this fails to converge to solution, logLogisticRegression samples lattice points throughout the parameter space.
#' It then uses the lattice point with minimal least-squares residual as an initial guess for the optimal parameters,
#' passes this guess to drm, and re-attempts the optimization. If this still fails, logLogisticRegression uses the
#' PatternSearch algorithm to fit a log-logistic curve to the data.
#'
#' @examples
#' dose <- c("0.0025","0.008","0.025","0.08","0.25","0.8","2.53","8") 
#' viability <- c("108.67","111","102.16","100.27","90","87","74","57")
#' computeAUC(dose, viability)
#' 
#' 
#' @param conc [vector] is a vector of drug concentrations.
#' @param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#' drug concentrations whose logarithms are in the corresponding entries of the log_conc, where viability 0
#' indicates that all cells died, and viability 1 indicates that the drug had no effect on the cells.   
#' @param density [vector] is a vector of length 3 whose components are the numbers of lattice points per unit
#' length along the HS-, E_inf-, and base-10 logarithm of the EC50-dimensions of the parameter space, respectively. 
#' @param step [vector] is a vector of length 3 whose entries are the initial step sizes in the HS, E_inf, and
#' base-10 logarithm of the EC50 dimensions, respectively, for the PatternSearch algorithm.
#' @param precision is a positive real number such that when the ratio of current step size to initial step
#' size falls below it, the PatternSearch algorithm terminates. A smaller value will cause LogisticPatternSearch
#' to take longer to complete optimization, but will produce a more accurate estimate for the fitted parameters.
#' @param lower_bounds [vector] is a vector of length 3 whose entries are the lower bounds on the HS, E_inf,
#' and base-10 logarithm of the EC50 parameters, respectively.
#' @param upper_bounds [vector] is a vector of length 3 whose entries are the upper bounds on the HS, E_inf,
#' and base-10 logarithm of the EC50 parameters, respectively.
#' @param scale is a positive real number specifying the shape parameter of the Cauchy distribution.
#' @param family [character], if "cauchy", uses MLE under an assumption of Cauchy-distributed errors
#' instead of sum-of-squared-residuals as the objective function for assessing goodness-of-fit of
#' dose-response curves to the data. Otherwise, if "normal", uses MLE with a gaussian assumption of errors
#' @param median_n If the viability points being fit were medians of measurements, they are expected to follow a median of \code{family}
#' distribution, which is in general quite different from the case of one measurement. Median_n is the number of measurements
#' the median was taken of. If the measurements are means of values, then both the Normal and the Cauchy distributions are stable, so means of 
#' Cauchy or Normal distributed variables are still Cauchy and normal respectively.
#' @param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data,
#' and that log10(EC50) should be returned instead of EC50.
#' @param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#' than a percentage, and that E_inf should be returned as a decimal rather than a percentage.
#' @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before 
#' curve-fitting is performed.
#' @param verbose [logical], if true, causes warnings thrown by the function to be printed.
#' @return A vector containing estimates for HS, E_inf, and EC50
#' @export
#' @importFrom stats optim dcauchy dnorm pcauchy rcauchy rnorm pnorm integrate

logLogisticRegression <- function(conc,
                                  viability,
                                  density = c(2, 10, 2),
                                  step = .5 / density,
                                  precision = 0.05,
                                  lower_bounds = c(0, 0, -6),
                                  upper_bounds = c(4, 1, 6),
                                  scale = 0.07,
                                  family = c("normal", "Cauchy"),
                                  median_n = 1,
                                  conc_as_log = FALSE,
                                  viability_as_pct = TRUE,
                                  trunc = TRUE,
                                  verbose = FALSE) {
  # guess <- .logLogisticRegressionRaw(conc, viability, density , step, precision, lower_bounds, upper_bounds, scale, Cauchy_flag, conc_as_log, viability_as_pct, trunc, verbose)


# .logLogisticRegressionRaw <- function(conc,
#                                   viability,
#                                   density = c(2, 10, 2),
#                                   step = .5 / density,
#                                   precision = 0.05,
#                                   lower_bounds = c(0, 0, -6),
#                                   upper_bounds = c(4, 1, 6),
#                                   scale = 0.07,
#                                   Cauchy_flag = FALSE,
#                                   conc_as_log = FALSE,
#                                   viability_as_pct = TRUE,
#                                   trunc = TRUE,
#                                   verbose = FALSE) {
  family <- match.arg(family)

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
  
  if (is.character(family) == FALSE) {
    print(family)
    stop("Cauchy flag is not a string.")
  }
  
  if (length(density) != 3){
    stop("Density parameter needs to have length of 3, for HS, Einf, EC50")
  }

  if (!median_n==as.integer(median_n)){
    stop("There can only be a integral number of samples to take a median of. Check your setting of median_n parameter, it is not an integer")
  }
  
  if (min(upper_bounds - lower_bounds) < 0) {
    print(rbind(lower_bounds, upper_bounds))
    stop("Upper bounds on parameters do not exceed lower bounds.")
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
  
  cleanData  <- sanitizeInput(conc=conc,
            viability=viability,
            conc_as_log = conc_as_log,
            viability_as_pct = viability_as_pct,
            trunc = trunc,
            verbose=verbose)

  log_conc <- cleanData[["log_conc"]]
  viability <- cleanData[["viability"]]
  
  
  #ATTEMPT TO REFINE GUESS WITH L-BFGS OPTIMIZATION
  # tryCatch(
  gritty_guess <- c(pmin(pmax(1, lower_bounds[1]), upper_bounds[1]),
                    pmin(pmax(min(viability), lower_bounds[2]), upper_bounds[2]),
                    pmin(pmax(log_conc[which.min(abs(viability - 1/2))], lower_bounds[3]), upper_bounds[3]))
  guess <- tryCatch(optim(par=gritty_guess,#par = sieve_guess,
                          fn = function(x) {.residual(log_conc,
                                                      viability,
                                                      pars = x,
                                                      n = median_n,
                                                      scale = scale,
                                                      family = family,
                                                      trunc = trunc)
                          },
                          lower = lower_bounds,
                          upper = upper_bounds,
                          method = "L-BFGS-B",
                          #control=list(maxit=100000)
                          ),#[[1]],
                    error = function(e) {
                      list("par"=gritty_guess, "convergence"=-1)
                    })
  # if(guess[["convergence"]]!=0)  print(guess[["message"]]); #print(results[["convergence"]])
  failed = guess[["convergence"]]!=0
  guess <- guess[["par"]]
  
  guess_residual <- .residual(log_conc,
                              viability,
                              pars = guess,
                              n = median_n,
                              scale = scale,
                              family = family,
                              trunc = trunc)
  

  #GENERATE INITIAL GUESS BY OBJECTIVE FUNCTION EVALUATION AT LATTICE POINTS

  gritty_guess_residual <- .residual(log_conc,
                                    viability,
                                    pars = gritty_guess,
                                    n = median_n,
                                    scale = scale,
                                    family = family,
                                    trunc = trunc)
  

  #CHECK SUCCESS OF L-BFGS OPTIMIZAITON AND RE-OPTIMIZE WITH A PATTERN SEARCH IF NECESSARY
  if (failed || any(is.na(guess)) || guess_residual >= gritty_guess_residual) {
    #GENERATE INITIAL GUESS BY OBJECTIVE FUNCTION EVALUATION AT LATTICE POINTS
    sieve_guess <- .meshEval(log_conc,
                           viability,
                           lower_bounds = lower_bounds,
                           upper_bounds = upper_bounds,
                           density = density,
                           n=median_n,
                           scale = scale,
                           family = family,
                           trunc = trunc)

    sieve_guess_residual <- .residual(log_conc,
                                    viability,
                                    pars = sieve_guess,
                                    n = median_n,
                                    scale = scale,
                                    family = family,
                                    trunc = trunc)

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
                                            n = median_n,
                                            scale = scale,
                                            family = family,
                                            trunc = trunc)
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
