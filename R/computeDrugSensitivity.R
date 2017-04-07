.calculateSensitivitiesStar <-
  function (pSets = list(), exps=NULL, cap=NA, na.rm=TRUE, area.type=c("Fitted","Actual"), nthread=1) {
    
    if (missing(area.type)) {
      area.type <- "Fitted"
    }
    if (is.null(exps)) {
      stop("expriments is empty!")
    }
    for (study in names(pSets)) {
      pSets[[study]]@sensitivity$profiles$auc_recomputed_star <- NA
    }
    if (!is.na(cap)) {
      trunc <- TRUE
    }else{
      trunc <- FALSE
    }
    
    for(i in 1:nrow(exps)) {
      ranges <- list()
      for (study in names(pSets)) {
        ranges[[study]] <- as.numeric(pSets[[study]]@sensitivity$raw[exps[i,study], ,"Dose"])
      }
      ranges <- .getCommonConcentrationRange(ranges)
      names(ranges) <- names(pSets)
      for(study in names(pSets)) {
      	myx <- as.numeric(pSets[[study]]@sensitivity$raw[exps[i, study],,"Dose"]) %in% ranges[[study]]
      	pSets[[study]]@sensitivity$raw[exps[i,study],!myx, ] <- NA
        
       }
    }
   	cl <- makeCluster(nthread)
    for(study in names(pSets)){

    	auc_recomputed_star <- unlist(parSapply(cl=cl, rownames(pSets[[study]]@sensitivity$raw), function(experiment, exps, study, dataset, area.type){
    		if(!experiment %in% exps[,study]){return(NA_real_)}
    		return(computeAUC(concentration=as.numeric(dataset[experiment,,1]), 
                        viability=as.numeric(dataset[experiment,,2]), 
 						trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE, area.type=area.type)/100)
 

    		}, exps = exps, study = study, dataset = pSets[[study]]@sensitivity$raw, area.type=area.type))
    	
    	pSets[[study]]@sensitivity$profiles$auc_recomputed_star <- auc_recomputed_star
    }
    stopCluster(cl)
    return(pSets)
  }

## This function computes AUC for the whole raw sensitivity data of a pset
.calculateFromRaw <- function(raw.sensitivity, cap=NA, nthread=1, family=c("normal", "Cauchy"), scale = 0.07, n = 1){
  family <- match.arg(family)

  AUC <- vector(length=dim(raw.sensitivity)[1])
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  
  IC50 <- vector(length=dim(raw.sensitivity)[1])
  names(IC50) <- dimnames(raw.sensitivity)[[1]]
  
  #pars <- logLogisticRegression(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"], conc_as_log=FALSE, viability_as_pct=TRUE, trunc=trunc)
  
  if (!is.na(cap)) {trunc <- TRUE}else{trunc <- FALSE}
  
  if (nthread ==1){
    pars <- lapply(names(AUC), function(exp, raw.sensitivity, family, scale, n) {
      if(length(grep("///", raw.sensitivity[exp, , "Dose"])) > 0 | all(is.na(raw.sensitivity[exp, , "Dose"]))) {
        NA
      } else{
        logLogisticRegression(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE, family=family, scale=scale, median_n=n)
        #computeAUC(concentration=raw.sensitivity[exp, , "Dose"], Hill_fit=Hill_fit, trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    },raw.sensitivity=raw.sensitivity, family = family, scale = scale, n = n)
    names(pars) <- dimnames(raw.sensitivity)[[1]]
    AUC <- unlist(lapply(names(pars), function(exp,raw.sensitivity, pars) {
      if(any(is.na(pars[[exp]]))) {
        NA
      } else{
        computeAUC(concentration=raw.sensitivity[exp, , "Dose"], Hill_fit=pars[[exp]], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    },raw.sensitivity=raw.sensitivity, pars=pars))
    IC50 <- unlist(lapply(names(pars), function(exp, pars) {
      if(any(is.na(pars[[exp]]))) {
        NA
      } else{
        computeIC50(Hill_fit=pars[[exp]], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    }, pars=pars))
  } else {
    pars <- parallel::mclapply(names(AUC), function(exp, raw.sensitivity, family, scale, n, trunc) {
      if(length(grep("///", raw.sensitivity[exp, , "Dose"])) > 0 | all(is.na(raw.sensitivity[exp, , "Dose"]))) {
        NA
      } else{
        logLogisticRegression(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE, family=family, scale=scale, median_n=n)
        #computeAUC(concentration=raw.sensitivity[exp, , "Dose"], Hill_fit=Hill_fit, trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    },raw.sensitivity=raw.sensitivity, family = family, scale = scale, n = n, trunc = trunc, mc.cores = nthread)
    names(pars) <- dimnames(raw.sensitivity)[[1]]
    AUC <- unlist(parallel::mclapply(names(pars), function(exp, raw.sensitivity, pars, trunc) {
      if(any(is.na(pars[[exp]]))) {
        NA
      } else{
        computeAUC(concentration=raw.sensitivity[exp, , "Dose"], Hill_fit=pars[[exp]], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    },raw.sensitivity=raw.sensitivity, pars=pars, trunc = trunc, mc.cores = nthread))
    IC50 <- unlist(parallel::mclapply(names(pars), function(exp, pars, trunc) {
      if(any(is.na(pars[[exp]]))) {
        NA
      } else{
        computeIC50(Hill_fit=pars[[exp]], trunc=trunc, conc_as_log=FALSE, viability_as_pct=TRUE)
      }
    }, pars=pars, trunc = trunc, mc.cores = nthread))
  }
  
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  names(IC50) <- dimnames(raw.sensitivity)[[1]]
  
  
  return(list("AUC"=AUC, "IC50"=IC50, "pars"=pars))
}


## This function computes intersected concentration range between a list of concentration ranges
.getCommonConcentrationRange <- function(doses)
{
  min.dose <- 0
  max.dose <- 10^100
  for(i in 1:length(doses))
  {
    min.dose <- max(min.dose, min(as.numeric(doses[[i]]), na.rm = TRUE), na.rm = TRUE)
    max.dose <- min(max.dose, max(as.numeric(doses[[i]]), na.rm = TRUE), na.rm = TRUE)
  }
  
  common.ranges <- list()
  for(i in 1:length(doses))
  {
    common.ranges[[i]] <- doses[[i]][
      which.min(abs(as.numeric(doses[[i]])-min.dose)):max(
        which(abs(as.numeric(doses[[i]]) - max.dose)==min(abs(as.numeric(doses[[i]]) - max.dose), na.rm=TRUE)))]
  }
  return(common.ranges)
}

## predict viability from concentration data and curve parameters
.Hill<-function(x, pars) {
  return(pars[2] + (1 - pars[2]) / (1 + (10 ^ x / 10 ^ pars[3]) ^ pars[1]))
}

## calculate residual of fit
.residual<-function(x, y, n, pars, scale = 0.07, family = c("normal", "Cauchy"), trunc = FALSE) {
  family <- match.arg(family)
  Cauchy_flag = (family == "Cauchy")
  if (Cauchy_flag == FALSE) {
    # return(sum((.Hill(x, pars) - y) ^ 2))
    diffs <- .Hill(x, pars)-y
    if (trunc == FALSE) {
      return(sum(-log(.dmednnormals(diffs, n, scale))))
    } else {
      down_truncated <- abs(y) >= 1
      up_truncated <- abs(y) <= 0
      
      # For up truncated, integrate the cauchy dist up until -diff because anything less gets truncated to 0, and thus the residual is -diff, and the prob
      # function becomes discrete
      # For down_truncated, 1-cdf(diffs) = cdf(-diffs) 
      
      return(sum(-log(.dmednnormals(diffs[!(down_truncated | up_truncated)], n, scale))) + sum(-log(.edmednnormals(-diffs[up_truncated | down_truncated], n, scale))))
      
    }
    
  } else {
    diffs <- .Hill(x, pars)-y
    if (trunc == FALSE) {
      # return(sum(-log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs / scale)) * (1 / 2 - 1 / pi * atan(diffs / scale)))))
      return(sum(-log(.dmedncauchys(diffs, n, scale))))
    } else {
      down_truncated <- abs(y) >= 1
      up_truncated <- abs(y) <= 0
      
      # For up truncated, integrate the cauchy dist up until -diff because anything less gets truncated to 0, and thus the residual is -diff, and the prob
      # function becomes discrete
      # For down_truncated, 1-cdf(diffs) = cdf(-diffs) 
      
      return(sum(-log(.dmedncauchys(diffs[!(down_truncated | up_truncated)], n, scale))) + sum(-log(.edmedncauchys(-diffs[up_truncated | down_truncated], n, scale))))
      
      # return(sum(log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))
      # * (1 / 2 - 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))),
      # -log(1 / 2 - 3 / (2 * pi) * atan((1 - diffs[down_truncated] - y[down_truncated]) / scale) + 2 / pi ^ 3 * (atan((1 - diffs[down_truncated] - y[down_truncated]) / scale)) ^ 3),
      # -log(-1 / 2 + 3 / (2 * pi) * atan((-diffs[up_truncated] - y[up_truncated]) / scale) - 2 / pi ^ 3 * (atan((- diffs[up_truncated] - y[up_truncated]) / scale)) ^ 3)))
    }
  }
}

## generate an initial guess for dose-response curve parameters by evaluating the residuals at different lattice points of the search space
.meshEval<-function(log_conc,
                    viability,
                    lower_bounds = c(0, 0, -6),
                    upper_bounds = c(4, 1, 6),
                    density = c(2, 10, 2),
                    scale = 0.07,
                    n = 1,
                    family = c("normal", "Cauchy"),
                    trunc = FALSE) {
  family <- match.arg(family)
  guess <- c(pmin(pmax(1, lower_bounds[1]), upper_bounds[1]),
                    pmin(pmax(min(viability), lower_bounds[2]), upper_bounds[2]),
                    pmin(pmax(log_conc[which.min(abs(viability - 1/2))], lower_bounds[3]), upper_bounds[3]))
  guess_residual<- .residual(log_conc,
                             viability,
                             pars = guess,
                             n=n,
                             scale = scale,
                             family = family,
                             trunc = trunc)
  for (i in seq(from = lower_bounds[1], to = upper_bounds[1], by = 1 / density[1])) {
    for (j in seq(from = lower_bounds[2], to = upper_bounds[2], by = 1 / density[2])) {
      for (k in seq(from = lower_bounds[3], to = upper_bounds[3], by = 1 / density[3])) {
        test_guess_residual <- .residual(log_conc,
                                         viability,
                                         pars = c(i, j, k),
                                         n=n,
                                         scale = scale,
                                         family = family,
                                         trunc = trunc)
        if(!is.finite(test_guess_residual)){
          warning(paste0(" Test Guess Residual is: ", test_guess_residual, "\n Other Pars: log_conc: ", paste(log_conc, collapse=", "), "\n Viability: ", paste(viability, collapse=", "), "\n Scale: ", scale, "\n Family: ", family, "\n Trunc ", trunc, "\n HS: ", i, ", Einf: ", j, ", logEC50: ", k, "\n n: ", n))
        }
        if(!length(test_guess_residual)){
          warning(paste0(" Test Guess Residual is: ", test_guess_residual, "\n Other Pars: log_conc: ", paste(log_conc, collapse=", "), "\n Viability: ", paste(viability, collapse=", "), "\n Scale: ", scale, "\n Family: ", family, "\n Trunc ", trunc, "\n HS: ", i, ", Einf: ", j, ", logEC50: ", k, "\n n: ", n))
        }
        if (test_guess_residual < guess_residual) {
          guess <- c(i, j, k)
          guess_residual <- test_guess_residual
        }
      }
    }
  }
  return(guess)
}

## get vector of interpolated concentrations for graphing purposes
.GetSupportVec <- function(x, output_length = 1001) {
  return(seq(from = min(x), to = max(x), length.out = output_length))
}
######## TODO ADD computationg from  being passed in params


#'  Fits dose-response curves to data given by the user
#'  and returns the AUC of the fitted curve, normalized to the length of the concentration range. 
#'
#'  @param concentration [vector] is a vector of drug concentrations.
#'
#'  @param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#'  drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
#'  of viability in the absence of any drug.
#'
#'  @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#'  curve-fitting is performed.
.computeAUCUnderFittedCurve <- function(concentration, viability, trunc=TRUE, verbose=FALSE) {
  
  # #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  # if (prod(is.finite(conc)) != 1) {
  #   print(conc)
  #   stop("Concentration vector contains elements which are not real numbers.")
  # }
  
  # if (prod(is.finite(viability)) != 1) {
  #   print(viability)
  #   stop("Viability vector contains elements which are not real numbers.")
  # }
  
  # if (is.logical(trunc) == FALSE) {
  #   print(trunc)
  #   stop("'trunc' is not a logical.")
  # }
  
  # if (length(conc) != length(viability)) {
  #   print(conc)
  #   print(viability)
  #   stop("Concentration vector is not of same length as viability vector.")
  # }
  
  # if (min(conc) < 0) {
  #   stop("Concentration vector contains negative data.")
  # }
  
  # if (min(viability) < 0 && verbose) {
  #   warning("Warning: Negative viability data.")
  # }
  
  # if (max(viability) > 100 && verbose) {
  #   warning("Warning: Viability data exceeds negative control.")
  # }
  
  # #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  # log_conc <- log10(conc)
  # viability <- viability / 100
  
  # if (trunc == TRUE) {
  #   viability[which(viability < 0)] <- 0
  #   viability[which(viability > 1)] <- 1
  # }
  log_conc <- concentration
  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(log_conc,
                                       viability,
                                       conc_as_log = TRUE,
                                       viability_as_pct = FALSE,
                                       trunc = trunc))
  x <- .GetSupportVec(log_conc)
  return(1 - trapz(x, .Hill(x, pars)) / (log_conc[length(log_conc)] - log_conc[1]))
}
#This function is being used in computeSlope 
.optimizeRegression <- function(x, y, x0 = -3, y0 = 100)
{
  beta1 = (sum(x * y) - y0 * sum(x)) / (sum(x * x) - x0 * sum(x))
  return(beta1)
}

updateMaxConc <- function(pSet){
  
  pSet@sensitivity$info$max.conc <- apply(pSet@sensitivity$raw[,,"Dose"], 1, max, na.rm=TRUE)
  return(pSet)
}

