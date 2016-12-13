




computeDSS <- function(concentration,
   viability,
   Hill_fit,
   t_param = 10,
   conc_as_log = FALSE,
   viability_as_pct = TRUE,
   trunc = TRUE,
   verbose = TRUE,
   dss_type=3
   #, ...
   ) {

  if(missing(concentration)){

    stop("The concentration values to integrate over must always be provided.")

}
if (missing(Hill_fit)) {

    Hill_fit <- logLogisticRegression(concentration,
      viability,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    cleanData <- sanitizeInput(conc=concentration, 
      Hill_fit=Hill_fit,
      conc_as_log = conc_as_log,
      viability_as_pct = viability_as_pct,
      trunc = trunc,
      verbose = verbose)
    pars <- cleanData[["Hill_fit"]]
    concentration <- cleanData[["log_conc"]]

} else {

  cleanData <- sanitizeInput(conc = concentration, 
    viability = viability,
    Hill_fit = Hill_fit,
    conc_as_log = conc_as_log,
    viability_as_pct = viability_as_pct,
    trunc = trunc,
    verbose = verbose)
  pars <- cleanData[["Hill_fit"]]
  concentration <- cleanData[["log_conc"]]
} 



# ic50 <- computeIC50(conc = concentration, Hill_fit = unlist(pars), conc_as_log = TRUE, viability_as_pct = FALSE)

# if (!is.finite(ic50) | ic50 > max(concentration)){
# 	return(0)
# }

if(pars[[3]] > max(concentration)){
	return(0)
}

if(viability_as_pct){
	t_param = t_param/100
}

x2 = max(concentration)
x1 = computeICn(conc = concentration, Hill_fit = unlist(pars), n = t_param, conc_as_log = TRUE, viability_as_pct = FALSE)
if(!is.finite(x1)){return(0)}

x1 <- max(x1, min(concentration))

AUC <- computeAUC(concentration = c(x1, x2), Hill_fit = unlist(pars), conc_as_log = TRUE, viability_as_pct = FALSE, verbose = verbose, trunc = trunc)


DSS <- ( AUC*(x2 - x1) - t_param * (x2 - x1))/((1 - t_param) * (max(concentration) - min(concentration))) * 100
if(dss_type == 1){
	return(DSS)
}
DSS <- DSS / (log10(1-pars[[2]])+2)
if(dss_type == 2){
	return(DSS)
}
DSS <- DSS * (x2 - x1)/(max(concentration) - min(concentration))
return(DSS)
}
