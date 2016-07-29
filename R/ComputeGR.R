GRRegression <- function(conc,
                  viability,
                  duration,
                  Hill_fit,
                  dbl_time,
                  conc_as_log = FALSE,
                  viability_as_pct = TRUE,
                  verbose = FALSE,
                  density = c(2, 10, 2),
                  step = .5 / density,
                  precision = 0.05,
                  lower_bounds = c(0, 0, -6),
                  upper_bounds = c(4, 1, 6),
                  scale = 0.07,
                  family = c("normal", "Cauchy"),
                  scale_01 = FALSE,
                  trunc=TRUE) { #If true, fits parameters to a transformed GR curve
  #with image [0, 1]. If false, fits parameters to Sorger's original [-1, 1]-image curve.
  
  #GRFit takes in dose-response data and a bunch of formatting parameters, then returns
  #the GHS, GEC_50, and G_inf values associated with them in accordance with the Sorger
  #paper. However, the definitions are corrected in accordance with my adjustment of the
  #relevant Sorger equations. While this may change the form of the equations, it does
  #not affect their intuitive meanings. G_inf is still the GR in the presence of arbitrarily
  #large drug concentration, GEC_50 is the dose that produces a half-minimal GR-value,
  #and GHS is the magnitude of the slope of the tangent to the log dose-response curve
  #when the drug concentration is GEC_50.
  
  #DO SANITY CHECKS ON INPUT
  # if(missing(concentration)){

  #   stop("The concentration values the drug was tested on must always be provided.")
  family <- match.arg(family)


if (missing(Hill_fit)) {

    Hill_fit <- logLogisticRegression(conc,
      viability,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose,
      density = density,
      step = step,
      precision = precision,
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      scale = scale,
      family = family
      )
    cleanData <- sanitizeInput(conc=conc, 
      Hill_fit=Hill_fit,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    Hill_fit <- cleanData[["Hill_fit"]]
    log_conc <- cleanData[["log_conc"]]
} else if (!missing(Hill_fit)){

  cleanData <- sanitizeInput(conc = conc, 
    viability = viability,
    Hill_fit = Hill_fit,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose)
  Hill_fit <- cleanData[["Hill_fit"]]
  # log_conc <- cleanData[["log_conc"]]
}
  
  if (missing(viability) && missing(Hill_fit)) {
    stop("Please enter viability data and/or Hill equation parameters.")
  }

  if(missing(duration)){
    stop("Cannot calculate GR without duration of experiment")
  }
  if(missing(dbl_time)){
    stop("Cannot calculate GR without cell doubling time")
  }
  
  tau <- duration / dbl_time
  
  #CALCULATE GR STATISTICS
  
  Ginf <- (Hill_fit[2]) ^ (1 / tau)
  GEC50 <- 10 ^ Hill_fit[3] * ((2 ^ tau - (1 + Ginf) ^ tau) / ((1 + Ginf) ^ tau - (2 * Ginf) ^ tau)) ^ (1 / Hill_fit[1])
  GHS <- (1 - Hill_fit[2]) / tau *
    (.Hill(log10(GEC50), Hill_fit)) ^ (1 / tau - 1) *
    1 / (1 + (GEC50 / 10 ^ Hill_fit[3]) ^ Hill_fit[1]) ^ 2 * Hill_fit[1] / 10 ^ Hill_fit[3] *
    (GEC50 / 10 ^ Hill_fit[3]) ^ (Hill_fit[1] - 1) * GEC50 * log(10)
  
  #CONVERT OUTPUT TO CONFORM TO FORMATTING PARAMETERS
  
  if (scale_01 == FALSE) {
    Ginf <- 2 * Ginf - 1
    GHS <- 2 * GHS
  }
  
  if (viability_as_pct == TRUE) {
    Ginf <- 100 * Ginf
    GHS <- 100 * GHS
  }
  
  if (conc_as_log == TRUE) {
    GEC50 <- log10(GEC50)
  }
  
  return(list(GHS = as.numeric(GHS), Ginf = as.numeric(Ginf), GEC50 = as.numeric(GEC50)))
  
}