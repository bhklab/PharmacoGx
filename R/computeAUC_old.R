# Return AUC (Area Under the drug response curve) for an experiment of a pSet by taking 
# its concentration and viability as input.
# 
# @param conc [vector] A concentration range that the AUC should be computed for that range.
# Concentration range by default considered as not logarithmic scaled.
# @param viability [vector] Viablities correspondant to the concentration range passed as first parameter.
# The range of viablity values by definition should be between 0 and 100. But the viabalities greater than
# 100 and lower than 0 are also accepted.
# @param trunc [binary] A flag that identify if the viabality values should be truncated to be in the
# range of (0,100)
# @param verbose [boolean] If 'TRUE' the function will retrun warnings and other infomrative messages.
# @import caTools
computeAUC_old <- function(conc, viability, 
                       conc_as_log = FALSE,
                       viability_as_pct = TRUE,
                       trunc=TRUE, 
                       verbose=TRUE, 
                       area.type=c("Fitted","Actual")) {
    cleanData <- sanitizeInput(conc, viability,
                             conc_as_log=conc_as_log,
                             viability_as_pct=viability_as_pct, 
                             trunc=trunc, verbose=verbose)
    log_conc <- cleanData[["log_conc"]]
    viability <- cleanData[["viability"]]

#   ii <- which(concentration == 0)
#   if(length(ii) > 0) {
#     concentration <- concentration[-ii]
#     viability <- viability[-ii]
#   }
  
  if(missing(area.type)){
    area.type <- "Fitted"
  }
  if(length(conc) < 2){
    return(NA)
  }
  if(area.type == "Actual"){
    # if(trunc) {viability = pmin(as.numeric(viability), 100); viability = pmax(as.numeric(viability), 0)}
    trapezoid.integral <- caTools::trapz(log10(as.numeric(conc) + 1) ,as.numeric(viability))
    AUC <- round(1- (trapezoid.integral/trapz(log10(as.numeric(conc)), rep(100, length(viability)))), digits=2)
  }else{
    AUC <- .computeAUCUnderFittedCurve(conc, viability, trunc)
  }
  return (AUC)
}
  
