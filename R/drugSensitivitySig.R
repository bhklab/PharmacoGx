################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:    
##      - data: gene expression data matrix
##            - drugpheno: sensititivity values fo thr drug of interest
##            - type: cell or tissue type for each experiment
##            - duration: experiment duration in hours
##      - batch: experiment batches
##            - single.type: Should the statitsics be computed for each cell/tissue type separately?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs:
## list of datafraes with the statistics for each gene, for each type
##
## Notes:    duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################


#' Creates a signature representing the association between gene expression (or
#' other molecular profile) and drug dose response, for use in drug sensitivity
#' analysis.
#' 
#' Given a Pharmacoset of the sensitivity experiment type, and a list of drugs, 
#' the function will compute a signature for the effect gene expression on the
#' molecular profile of a cell. The function returns the estimated coefficient, 
#' the t-stat, the p-value and the false discovery rate associated with that 
#' coefficient, in a 3 dimensional array, with genes in the first direction, 
#' drugs in the second, and the selected return values in the third.
#' 
#' @examples
#' data(GDSCsmall)
#' drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", 
#'              nthread=1, features = featureNames(GDSCsmall, "rna")[1])
#' print(drug.sensitivity)
#' 
#' @param pSet [PharmacoSet] a PharmacoSet of the perturbation experiment type
#' @param mDataType [character] which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv
#' @param drugs [character] a vector of drug names for which to compute the
#'   signatures. Should match the names used in the PharmacoSet.
#' @param features [character] a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in PharmacoSet.
#' @param nthread [numeric] if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues [character] Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene drug pair?
#' @param sensitivity.measure [character] which measure of the drug dose 
#'   sensitivity should the function use for its computations? The current
#'   choices are 'ic50_published', 'auc_published', 'ic50_recomputed',
#'   'auc_recomputed'.
#' @param molecular.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line molecular profile measurements? 
#' @param sensitivity.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line sensitivity measurements? 
#' @param sensitivity.cutoff Allows to provide upper and lower bounds to
#'   sensitivity measures in the cases where the values exceed physical values
#'   due to numerical or other errors.
#' @param verbose [boolean] 'TRUE' if the warnings and other infomrative message shoud be displayed
#' @return [list] a 3D array with genes in the first dimension, drugs in the
#'   second, and return values in the third.
#' @export
#' @import parallel

drugSensitivitySig <- function(pSet, mDataType, drugs, features, sensitivity.measure=c("ic50_published", "auc_published", "ic50_recomputed", "auc_recomputed", "auc_recomputed_star"), molecular.summary.stat=c("mean", "median", "first", "last", "or", "and"), sensitivity.summary.stat=c("mean", "median", "first", "last"), returnValues=c("estimate", "pvalue", "fdr"), sensitivity.cutoff, nthread=1, verbose=TRUE) {
    
    ### This function needs to: Get a table of AUC values per cell line / drug
    ### Be able to recompute those values on the fly from raw data if needed to change concentration
    ### Be able to choose different summary methods on fly if needed (need to add annotation to table to tell what summary method previously used)
    ### Be able to extract genomic data 
    ### Run rankGeneDrugSens in parallel at the drug level
    ### Return matrix as we had before
  
  sensitivity.measure <- match.arg(sensitivity.measure)
  molecular.summary.stat <- match.arg(molecular.summary.stat)
  sensitivity.summary.stat <- match.arg(sensitivity.summary.stat)
  
  
  if (!(sensitivity.measure %in% colnames(sensitivityProfiles(pSet)))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", pSet@annotation$name, paste(colnames(sensitivityProfiles(pSet)), collapse=", ")))
  }
  
  if (!(mDataType %in% names(pSet@molecularProfiles))) {
    stop (sprintf("Invalid mDataType for %s, choose among: %s", pSet@annotation$name, paste(names(pSet@molecularProfiles), collapse=", ")))
  }
  switch (Biobase::annotation(pSet@molecularProfiles[[mDataType]]),
    "mutation" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop ("Molecular summary statistic for mutation must be either 'or' or 'and'")
      }
    },
    "fusion" = {
      if (!is.element(molecular.summary.stat, c("or", "and"))) {
        stop ("Molecular summary statistic for fusion must be either 'or' or 'and'")
      }
    },
    "rna" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    "cnv" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for cnv must be either 'mean', 'median', 'first' or 'last'")
      }
    },
    stop (sprintf("No summary statistic for %s has been implemented yet", Biobase::annotation(pSet@molecularProfiles[[mDataType]])))
  )
  
  if (!is.element(sensitivity.summary.stat, c("mean", "median", "first", "last"))) {
    stop ("Sensitivity summary statistic for sensitivity must be either 'mean', 'median', 'first' or 'last'")
  }
  
  if (missing(sensitivity.cutoff)) {
    sensitivity.cutoff <- NA
  }
  if (missing(drugs)){
    drugn <- drugNames(pSet)
  } else {
    drugn <- drugs
  }
  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  
  if (missing(features)) {
    features <- rownames(featureInfo(pSet, mDataType))
  } else {
    fix <- is.element(features, rownames(featureInfo(pSet, mDataType)))
    if (verbose && !all(fix)) {
      warning (sprintf("%i/%i features can be found", sum(fix), length(features)))
    }
    features <- features[fix]
  }
  
    drugpheno.all <- t(summarizeSensitivityProfiles(pSet, sensitivity.measure=sensitivity.measure, summary.stat=sensitivity.summary.stat, verbose=verbose))
    dix <- is.element(drugn, colnames(drugpheno.all))
  if (verbose && !all(dix)) {
    warning (sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
  }
  if (!any(dix)) {
    stop("None of the drugs were found in the dataset")
  }
  drugn <- drugn[dix]
  
    pSet@molecularProfiles[[mDataType]] <- summarizeMolecularProfiles(pSet=pSet, mDataType=mDataType, summary.stat=molecular.summary.stat, verbose=verbose)
    drugpheno.all <- drugpheno.all[phenoInfo(pSet, mDataType)[ ,"cellid"], , drop=FALSE]
    
    type <- as.factor(cellInfo(pSet)[phenoInfo(pSet, mDataType)[ ,"cellid"], "tissueid"]) 
  batch <- phenoInfo(pSet, mDataType)[, "batchid"]
  batch[!is.na(batch) & batch == "NA"] <- NA
  batch <- as.factor(batch)
    names(batch) <- phenoInfo(pSet, mDataType)[ , "cellid"]
    batch <- batch[rownames(drugpheno.all)]
    # duration <- sensitivityInfo(pSet)[,"duration_h"]
  ## compute drug sensitivity signatures
  if (verbose) {
    message("Computation of drug sensitivity signatures...")
  }
  splitix <- parallel::splitIndices(nx=length(drugn), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <-  parallel::mclapply(splitix, function(x, drugn, expr, drugpheno, type, batch) {
  res <- NULL
  for(i in drugn[x]) {
  ## using a linear model (x ~ concentration + cell + batch)
      dd <- drugpheno[ , i]
      if(!is.na(sensitivity.cutoff)) {
            dd <- factor(ifelse(dd > sensitivity.cutoff, 1, 0), levels=c(0, 1))
      }
      rr <- rankGeneDrugSensitivity(data=expr, drugpheno=dd, type=type, batch=batch, single.type=FALSE, nthread=1, verbose=verbose)
      res <- c(res, list(rr$all))
    }
    names(res) <- drugn[x]
    return(res)
  }, drugn=drugn, expr=t(molecularProfiles(pSet, mDataType)[features, , drop=FALSE]), drugpheno=drugpheno.all, type=type, batch=batch)
  res <- do.call(c, mcres)
  res <- res[!sapply(res, is.null)]
  drug.sensitivity <- array(NA, dim=c(nrow(featureInfo(pSet, mDataType)[features,, drop=FALSE]), length(res), ncol(res[[1]])), dimnames=list(rownames(featureInfo(pSet, mDataType)[features,]), names(res), colnames(res[[1]])))
  for(j in 1:ncol(res[[1]])) {
    ttt <- sapply(res, function(x, j, k) {
      xx <- array(NA, dim=length(k), dimnames=list(k))
      xx[rownames(x)] <- x[ , j, drop=FALSE]
      return (xx)
    }, j=j, k=rownames(featureInfo(pSet, mDataType)[features,, drop=FALSE]))
    drug.sensitivity[rownames(featureInfo(pSet, mDataType)[features,, drop=FALSE]), names(res), j] <- ttt
  }

drug.sensitivity <- PharmacoSig(drug.sensitivity, PSetName = pSetName(pSet), Call ="as.character(match.call())", SigType='Sensitivity')

    return(drug.sensitivity)
}
