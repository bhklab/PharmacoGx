#' @describeIn computeICn Returns the IC50 of a Drug Dose response curve
#' @export
computeIC50 <- function(concentration,
                       viability,
                       Hill_fit,
                       conc_as_log = FALSE,
                       viability_as_pct = TRUE, 
                       verbose = TRUE,
                       trunc = TRUE) {

  return(computeICn(concentration = concentration, 
                    viability = viability, 
                    Hill_fit = Hill_fit,
                    n = ifelse(viability_as_pct, 50, .5),
                    conc_as_log = conc_as_log,
                    viability_as_pct = viability_as_pct, 
                    verbose=TRUE,
                    trunc=TRUE))
}
