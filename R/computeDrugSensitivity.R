########################
## Mark Freeman, Zhaleh Safikhani, Petr Smirnov, Benjamin Haibe-Kains
## All rights Reserved
## May 26, 2015
########################


.calculateSensitivities <-
  function(pSets = list(), cellMatch=NULL, drugMatch=NULL, cap=NA, na.rm=TRUE){
#  require(caTools) || stop("Library caTools is not available!") # trapezoid function 
  
  cell.match.flag <- ifelse(is.null(cellMatch), T, F)
  drug.match.flag <- ifelse(is.null(drugMatch), T, F)
  
  if(cell.match.flag){cellMatch <- pSets[[1]]@curation$cell$unique.cellid}
  if(drug.match.flag){drugMatch <- pSets[[1]]@curation$drug$unique.drugid}
  
  for (i in 1:length(pSets))
  {
    sensitivityData(pSets[[i]]) <- cbind(sensitivityData(pSets[[i]]), "auc_recomputed_star" = NA)
    if(cell.match.flag){cellMatch <- intersect(cellMatch, pSets[[i]]@curation$cell$unique.cellid)}
    if(drug.match.flag){drugMatch <- intersect(drugMatch, pSets[[i]]@curation$drug$unique.drugid)}  
  }
  if (!is.na(cap)) {trunc <- TRUE}else{trunc <- FALSE}
  
  for(x in 1:length(cellMatch))
  {
      for(d in 1: length(drugMatch))
      {
        ranges <- list()
        exp_id <- list()
        flag <- F
        for(i in 1:length(pSets))
        {
          study <- names(pSets)[i]
          #exp_id[[i]] <- paste(drugMatch[d,study], cellMatch[x,paste0(study,".cellid")], sep = "_")
          exp_id[[i]] <- paste(x, d, sep = "_")
          if(is.null(exp_id[[i]]) | length(exp_id[[i]]) == 0){ stop("Studies should have names equal to the ones used in matching data frames")}
          if(exp_id[[i]] %in% dimnames(pSets[[i]]@sensitivity$raw)[[1]])
          {
            ranges[[i]] <- as.numeric(pSets[[i]]@sensitivity$raw[exp_id[[i]],,"Dose"])
          }else{
            flag <- T
          }
        }
        if(!flag)
        {
          ranges <- .getCommonConcentrationRange(ranges)
          for(i in 1:length(pSets))
          {
            sensitivityData(pSets[[i]])[exp_id[[i]], "auc_recomputed_star"] <- computeAUC(na.omit(ranges[[i]]), na.omit(pSets[[i]]@sensitivity$raw[exp_id[[i]],which(pSets[[i]]@sensitivity$raw[exp_id[[i]],,"Dose"] %in% ranges[[i]]),"Viability"]), trunc)
          }
        }
      } 
  }
  
  #pSet <- updateTables(pSet,summary=pSet@table.summary["auc_recomputed"], measures="auc_recomputed", verbose=FALSE)
  return(pSets)
}

## This function computes sensitivity for the whole raw sensitivity data of a pset
.calculateFromRaw <- function(raw.sensitivity, dose.range=vector(), cap=NA, na.rm=TRUE, nbcore = 1){
  #require(caTools) || stop("Library magicaxis is not available!")
  
  
  AUC <- vector(length=dim(raw.sensitivity)[1])
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  
  
  if (!is.na(cap)) {trunc <- TRUE}else{trunc <- FALSE}
  
  if(na.rm){
    if (nbcore ==1){
      AUC <- unlist(lapply(names(AUC), function(exp){computeAUC(na.omit(raw.sensitivity[exp,,"Dose"]), na.omit(raw.sensitivity[exp,,"Viability"]), trunc)}))
    }else{ 
      AUC <- unlist(parallel::mclapply(names(AUC), function(exp){computeAUC(na.omit(raw.sensitivity[exp,,"Dose"]), na.omit(raw.sensitivity[exp,,"Viability"]), trunc)}, mc.cores = nbcore))
    }
  }else{
    if (nbcore ==1){
      AUC <- unlist(lapply(names(AUC), function(exp){computeAUC(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}))
    }else{            
      AUC <- unlist(parallel::mclapply(names(AUC)[1:10], function(exp){computeAUC(raw.sensitivity[exp,,"Dose"], raw.sensitivity[exp,,"Viability"], trunc)}, mc.cores = nbcore))
    }
  }
  
  names(AUC) <- dimnames(raw.sensitivity)[[1]]
  
  return(list("AUC" = AUC))
}


## This function computes intersected concentration range between a list of concentration ranges
.getCommonConcentrationRange <- function(doses)
{
  min.dose <- 0
  max.dose <- 10^100
  for(i in 1:length(doses))
  {
    min.dose <- max(min.dose, min(as.numeric(doses[[i]]), na.rm = T), na.rm = T)
    max.dose <- min(max.dose, max(as.numeric(doses[[i]]), na.rm = T), na.rm = T)
  }
  
  common.ranges <- list()
  for(i in 1:length(doses))
  {
    common.ranges[[i]] <- doses[[i]][which.min(abs(as.numeric(doses[[i]])-min.dose)):which.min(abs(as.numeric(doses[[i]])-max.dose))]
  }
  return(common.ranges)
}
#HIDDEN FUNCTIONS

#PREDICT VIABILITY FROM CONCENTRATION DATA AND CURVE PARAMETERS
.Hill<-function(x, pars) {
    return(pars[2] + (1 - pars[2]) / (1 + (10 ^ x / 10 ^ pars[3]) ^ pars[1]))
}

#CALCULATE RESIDUAL OF FIT
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

#GENERATE AN INITIAL GUESS FOR DOSE-RESPONSE CURVE PARAMETERS BY EVALUATING THE RESIDUALS AT DIFFERENT LATTICE POINTS OF THE SEARCH SPACE
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

#GET VECTOR OF INTERPOLATED CONCENTRATIONS FOR GRAPHING PURPOSES
.GetSupportVec <- function(x, output_length = 1001) {
  return(seq(from = min(x), to = max(x), length.out = output_length))
}
