#' @importFrom stats na.omit
.calculateSensitivitiesStar <-
  function (pSets = list(), exps=NULL, cap=NA, na.rm=TRUE, area.type=c("Fitted","Actual")) {
  
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
        pSets[[study]]@sensitivity$profiles[exps[i,study], "auc_recomputed_star"] <- computeAUC(concentration=as.numeric(na.omit(ranges[[study]])), 
                                                                                 viability=as.numeric(na.omit(pSets[[study]]@sensitivity$raw[exps[i, study],which(as.numeric(pSets[[study]]@sensitivity$raw[exps[i, study],,"Dose"]) %in% ranges[[study]]),"Viability"])), 
                                                                                 trunc, area.type)
      }
  }
  return(pSets)
}

## This function computes AUC for the whole raw sensitivity data of a pset
.calculateFromRaw <- function(raw.sensitivity, cap=NA, nthread=1){
  
  AUC <- vector(length=dim(raw.sensitivity)[1])
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  
  IC50 <- vector(length=dim(raw.sensitivity)[1])
  names(IC50) <- dimnames(raw.sensitivity)[[1]]
  
  
  
  if (!is.na(cap)) {trunc <- TRUE}else{trunc <- FALSE}
  
   if (nthread ==1){
    AUC <- unlist(lapply(names(AUC), function(exp){computeAUC(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}))
    IC50 <- unlist(lapply(names(IC50), function(exp){computeIC50(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}))
    
  }else{ 
    AUC <- unlist(parallel::mclapply(names(AUC), function(exp){computeAUC(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}, mc.cores = nthread))
    IC50 <- unlist(parallel::mclapply(names(IC50), function(exp){computeIC50(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}, mc.cores = nthread))
  }

  
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  names(IC50) <- dimnames(raw.sensitivity)[[1]]
  
  
  return(list("AUC"=AUC, "IC50"=IC50))
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
    common.ranges[[i]] <- doses[[i]][which.min(abs(as.numeric(doses[[i]])-min.dose)):which.min(abs(as.numeric(doses[[i]])-max.dose))]
  }
  return(common.ranges)
}

## predict viability from concentration data and curve parameters
.Hill<-function(x, pars) {
    return(pars[2] + (1 - pars[2]) / (1 + (10 ^ x / 10 ^ pars[3]) ^ pars[1]))
}

## calculate residual of fit
.residual<-function(x, y, pars, scale = 0.07, Cauchy_flag = TRUE, trunc = FALSE) {
    if (Cauchy_flag == FALSE) {
        return(sum((.Hill(x, pars) - y) ^ 2))
    } else {
        diffs <- .Hill(x, pars)-y
        if (trunc == FALSE) {
            return(sum(-log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs / scale)) * (1 / 2 - 1 / pi * atan(diffs / scale)))))
        } else {
            down_truncated <- which(abs(y) > 1)
            up_truncated <- which(abs(y) < 0)
            return(sum(log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))
            * (1 / 2 - 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))),
            -log(1 / 2 - 3 / (2 * pi) * atan((1 - diffs[down_truncated] - y[down_truncated]) / scale) + 2 / pi ^ 3 * (atan((1 - diffs[down_truncated] - y[down_truncated]) / scale)) ^ 3),
            -log(-1 / 2 + 3 / (2 * pi) * atan((-diffs[up_truncated] - y[up_truncated]) / scale) - 2 / pi ^ 3 * (atan((- diffs[up_truncated] - y[up_truncated]) / scale)) ^ 3)))
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
Cauchy_flag = TRUE,
trunc = FALSE) {
    guess_residual<-Inf
    for (i in seq(from = lower_bounds[1], to = upper_bounds[1], by = 1 / density[1])) {
        for (j in seq(from = lower_bounds[2], to = upper_bounds[2], by = 1 / density[2])) {
            for (k in seq(from = lower_bounds[3], to = upper_bounds[3], by = 1 / density[3])) {
                test_guess_residual <- .residual(log_conc,
                viability,
                pars = c(i, j, k),
                scale = scale,
                Cauchy_flag = Cauchy_flag,
                trunc = trunc)
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
#'  Fits dose-response curves to data given by the user
#'  and returns the AUC of the fitted curve, normalized to the length of the concentration range.
#'
#'  @param conc [vector] is a vector of drug concentrations.
#'
#'  @param viability [vector] is a vector whose entries are the viability values observed in the presence of the
#'  drug concentrations whose logarithms are in the corresponding entries of the log_conc, expressed as percentages
#'  of viability in the absence of any drug.
#'
#'  @param trunc [logical], if true, causes viability data to be truncated to lie between 0 and 1 before
#'  curve-fitting is performed.
.computeAUCUnderFittedCurve <- function(conc, viability, trunc=TRUE, verbose=FALSE) {
  
  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
  if (prod(is.finite(conc)) != 1) {
    print(conc)
    stop("Concentration vector contains elements which are not real numbers.")
  }
  
  if (prod(is.finite(viability)) != 1) {
    print(viability)
    stop("Viability vector contains elements which are not real numbers.")
  }
  
  if (is.logical(trunc) == FALSE) {
    print(trunc)
    stop("'trunc' is not a logical.")
  }
  
  if (length(conc) != length(viability)) {
    print(conc)
    print(viability)
    stop("Concentration vector is not of same length as viability vector.")
  }
  
  if (min(conc) < 0) {
    stop("Concentration vector contains negative data.")
  }
  
  if (min(viability) < 0 && verbose) {
    warning("Warning: Negative viability data.")
  }
  
  if (max(viability) > 100 && verbose) {
    warning("Warning: Viability data exceeds negative control.")
  }
  
  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
  log_conc <- log10(conc)
  viability <- viability / 100
  
  if (trunc == TRUE) {
    viability[which(viability < 0)] <- 0
    viability[which(viability > 1)] <- 1
  }
  
  #FIT CURVE AND CALCULATE IC50
  pars <- unlist(logLogisticRegression(log_conc,
                                       viability,
                                       conc_as_log = TRUE,
                                       viability_as_pct = FALSE,
                                       trunc = trunc))
  x <- .GetSupportVec(log_conc)
  return(1 - trapz(x, .Hill(x, pars)) / (log_conc[length(log_conc)] - log_conc[1]))
}

