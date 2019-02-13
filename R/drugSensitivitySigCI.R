################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:    
##      - data: gene expression data matrix
##            - drugpheno: sensitivity values for the drug of interest
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
#'              nthread=1, features = fNames(GDSCsmall, "rna")[1])
#' print(drug.sensitivity)
#' 
#' @param pSet [PharmacoSet] a PharmacoSet of the perturbation experiment type
#' @param mDataType [character] which one of the molecular data types to use
#'   in the analysis, out of dna, rna, rnaseq, snp, cnv
#' @param drugs [character] a vector of drug names for which to compute the
#'   signatures. Should match the names used in the PharmacoSet.
#' @param features [character] a vector of features for which to compute the
#'   signatures. Should match the names used in correspondant molecular data in PharmacoSet.
#' @param cells [character] allows choosing exactly which cell lines to include for the signature fitting. 
#'   Should be a subset of cellNames(pSet)
#' @param tissues [character] a vector of which tissue types to include in the signature fitting. 
#'   Should be a subset of cellInfo(pSet)$tissueid
#' @param nthread [numeric] if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues [character] Which of estimate, t-stat, p-value and fdr
#'   should the function return for each gene drug pair?
#' @param sensitivity.measure [character] which measure of the drug dose 
#'   sensitivity should the function use for its computations? Use the 
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param molecular.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line molecular profile measurements? 
#' @param sensitivity.summary.stat What summary statistic should be used to
#'   summarize duplicates for cell line sensitivity measurements? 
#' @param sensitivity.cutoff [numeric] Allows the user to binarize the sensitivity data using this threshold.
#' @param standardize [character] One of "SD", "rescale", or "none", for the form of standardization of
#'   the data to use. If "SD", the the data is scaled so that SD = 1. If rescale, then the data is scaled so that the 95%
#'   interquantile range lies in [0,1]. If none no rescaling is done. 
#' @param molecular.cutoff Allows the user to binarize the sensitivity data using this threshold. 
#' @param molecular.cutoff.direction [character] One of "less" or "greater", allows to set direction of binarization. 
#' @param verbose [boolean] 'TRUE' if the warnings and other infomrative message shoud be displayed
#' @param ... additional arguments not currently fully supported by the function  
#' @return [list] a 3D array with genes in the first dimension, drugs in the
#'   second, and return values in the third.
#' @export
#' @import parallel
#' @importFrom fastCI justFastCI combineCI.2
drugSensitivitySigCI <- function(pSet,
 mDataType,
 drugs,
 features,
 cells, 
 tissues,
 sensitivity.measure = "auc_recomputed", 
 molecular.summary.stat = c("mean", "median", "first", "last", "or", "and"), 
 sensitivity.summary.stat = c("mean", "median", "first", "last"), 
 returnValues = c("estimate", "pvalue", "fdr"),
 sensitivity.cutoff, standardize = c("rescale", "SD", "none"),
 molecular.cutoff = NA,
 molecular.cutoff.direction = c("less", "greater"),
 nthread = 1,
 verbose=TRUE,
 ...) {
  
  ### This function needs to: Get a table of AUC values per cell line / drug
  ### Be able to recompute those values on the fly from raw data if needed to change concentration
  ### Be able to choose different summary methods on fly if needed (need to add annotation to table to tell what summary method previously used)
  ### Be able to extract genomic data 
  ### Run rankGeneDrugSens in parallel at the drug level
  ### Return matrix as we had before
  
  #sensitivity.measure <- match.arg(sensitivity.measure)
  molecular.summary.stat <- match.arg(molecular.summary.stat)
  sensitivity.summary.stat <- match.arg(sensitivity.summary.stat)
  standardize <- match.arg(standardize)
  molecular.cutoff.direction <- match.arg(molecular.cutoff.direction)
  dots <- list(...)
  ndots <- length(dots)
  
  
  if (!all(sensitivity.measure %in% colnames(sensitivityProfiles(pSet)))) {
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
    "rnaseq" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
    }},
    "isoform" = {
      if (!is.element(molecular.summary.stat, c("mean", "median", "first", "last"))) {
        stop ("Molecular summary statistic for rna must be either 'mean', 'median', 'first' or 'last'")
    }},
    stop (sprintf("No summary statistic for %s has been implemented yet", Biobase::annotation(pSet@molecularProfiles[[mDataType]])))
  )
  
  if (!is.element(sensitivity.summary.stat, c("mean", "median", "first", "last"))) {
    stop ("Sensitivity summary statistic for sensitivity must be either 'mean', 'median', 'first' or 'last'")
  }
  
  if (missing(sensitivity.cutoff)) {
    sensitivity.cutoff <- NA
  }
  if (missing(drugs)){
    drugn <- drugs <- drugNames(pSet)
  } else {
    drugn <- drugs
  }

  if (missing(cells)){
    celln <- cells <- cellNames(pSet)
  } else {
    celln <- cells
  }

  availcore <- parallel::detectCores()
  if ( nthread > availcore) {
    nthread <- availcore
  }
  

  
  if(is.null(dots[["sProfiles"]])){
    drugpheno.all <- lapply(sensitivity.measure, function(sensitivity.measure) {
      
      return(t(summarizeSensitivityProfiles(pSet,
        sensitivity.measure = sensitivity.measure,
        summary.stat = sensitivity.summary.stat,
        verbose = verbose)))
      
    })} else {
      sProfiles <- dots[["sProfiles"]]
      drugpheno.all <- list(t(sProfiles))
    }
    
    dix <- is.element(drugn, do.call(colnames, drugpheno.all))
    if (verbose && !all(dix)) {
      warning (sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
    }
    if (!any(dix)) {
      stop("None of the drugs were found in the dataset")
    }
    drugn <- drugn[dix]

    cix <- is.element(celln, do.call(rownames, drugpheno.all))
    if (verbose && !all(cix)) {
      warning (sprintf("%i/%i cells can be found", sum(cix), length(celln)))
    }
    if (!any(cix)) {
      stop("None of the cells were found in the dataset")
    }
    celln <- celln[cix]
    
    if(!missing(tissues)){
      celln <- celln[cellInfo(pSet)[celln,"tissueid"] %in% tissues]
    } else {
      tissues <- unique(cellInfo(pSet)$tissueid)
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

    if(!is.null(dots[["mProfiles"]])){
      mProfiles <- dots[["mProfiles"]]
      Biobase::exprs(pSet@molecularProfiles[[mDataType]]) <- mProfiles[features, colnames(pSet@molecularProfiles[[mDataType]]), drop = FALSE]
      
    }
    pSet@molecularProfiles[[mDataType]] <- summarizeMolecularProfiles(pSet = pSet,
      mDataType = mDataType,
      summary.stat = molecular.summary.stat,
      binarize.threshold = molecular.cutoff,
      binarize.direction = molecular.cutoff.direction,
      verbose = verbose)[features, ]

    drugpheno.all <- lapply(drugpheno.all, function(x) {x[intersect(phenoInfo(pSet, mDataType)[ ,"cellid"], celln), , drop = FALSE]})
    
    molcellx <- phenoInfo(pSet, mDataType)[ ,"cellid"] %in% celln

    type <- as.factor(cellInfo(pSet)[phenoInfo(pSet, mDataType)[molcellx,"cellid"], "tissueid"]) 
    # batch <- phenoInfo(pSet, mDataType)[molcellx, "batchid"]
    # batch[!is.na(batch) & batch == "NA"] <- NA
    # batch <- as.factor(batch)
    # names(batch) <- phenoInfo(pSet, mDataType)[molcellx, "cellid"]
    # batch <- batch[rownames(drugpheno.all[[1]])]
    data <- t(Biobase::exprs(pSet@molecularProfiles[[mDataType]])[,molcellx])

    ## split task up by: gene, drug and tissue

    drugpheno <- drugpheno.all[[1]]
    
    if(verbose){
      message("Calculating Tissue Specific Statistics")
    }
    # browser()
    resLists <- foreach(dn = drugn) %:% 
      foreach(tissue = unique(type)) %dopar% {
        dp <- drugpheno[type == tissue, dn]
        dt <- data[type == tissue,]
        rankGeneDrugSensitivityCI(dt, dp, verbose=FALSE) 
      }

    #
    if(verbose){
      message("Computing Pan Cancer results")
    }

    resListsPanCancer <- foreach(dn = drugn) %dopar% {
        dp <- drugpheno[, dn]
        dt <- data
        rankGeneDrugSensitivityCI(dt, dp, verbose=FALSE) 
      }


    resLists <- lapply(resLists, function(x) return(abind(x, along = 3)))
    res <- abind(resLists, along = -1)
    dimnames(res)[[1]] <- drugn
    dimnames(res)[[4]] <- unique(type)
    
    pancancer.res <- abind(resListsPanCancer, along = -1)

    dimnames(pancancer.res)[[1]] <- drugn

    if(verbose){
      message("Assessing Significance")
    }

    # Switches to normal approx in case of 100 or more samples
    maxN <- max(pmin(as.vector(apply(res[,2,,], c(1,2),sum, na.rm=TRUE)), 100, na.rm=TRUE))

    nullTable <- makeTableUpToN(maxN)

    ## First calculate pan-cancer

    pancancer.p <- apply(pancancer.res, c(1,3), function(x){
      if(!is.finite(x[2])){
        return(NA_real_)
      }
      if(x[2] < 100){
        return(getCIPvals(cumsum(nullTable[[x[2]]]), x[1]*choose(x[2],2)))
      } else {
        pars.out <- computeExpectedApproximation(x[2])
        conc.out <- x[1]*choose(x[2],2)
        if(conc.out < pars.out[1]){
          prob <- pnorm(conc.out + 1,mean = pars.out[1], sd = pars.out[2])
          p.out <- 2*prob
        } else {
          prob <- pnorm(conc.out - 1,mean = pars.out[1], sd = pars.out[2], lower.tail = FALSE)
          p.out <- 2*prob
        }
        return(p.out)
      }
    })
    adjusted.pancancer.p <- p.adjust(pancancer.p, method = "fdr")
    dim(adjusted.pancancer.p) <- dim(pancancer.p)
    pancancer.res <- abind(pancancer.res, "p" = pancancer.p, "fdr" = adjusted.pancancer.p, along=2)


    ## Calculating p-values for combined tests
    tissue.stratified <- apply(res, c(1,3), combineCI.2, nullTable = nullTable)
    adjusted.p <- p.adjust(tissue.stratified[2,,], method = "fdr")
    dim(adjusted.p) <- dim(res)[c(1,3)]
    tissue.stratified <- abind(tissue.stratified, "fdr" = adjusted.p, along = 1)
    tissue.stratified <- aperm(tissue.stratified, c(2,1,3))
    # Tissue Specific P Values
    tissue.p <- apply(res, c(1,3,4), function(x){
      if(!is.finite(x[2]) | x[2] == 0){
        return(NA_real_)
      }
      if(x[2] < 100){
        return(getCIPvals(cumsum(nullTable[[x[2]]]), x[1]*choose(x[2],2)))
      } else {
        pars.out <- computeExpectedApproximation(x[2])
        conc.out <- x[1]*choose(x[2],2)
        if(conc.out < pars.out[1]){
          prob <- pnorm(conc.out + 1,mean = pars.out[1], sd = pars.out[2])
          p.out <- 2*prob
        } else {
          prob <- pnorm(conc.out - 1,mean = pars.out[1], sd = pars.out[2], lower.tail = FALSE)
          p.out <- 2*prob
        }
        return(p.out)
      }
    })
    adjusted.tissue.p <- p.adjust(tissue.p, method = "fdr")
    dim(adjusted.tissue.p) <- dim(tissue.p)
    res <- abind(res, "p" = tissue.p, "fdr" = adjusted.tissue.p, along=2)

    return(list("Pan Cancer" = pancancer.res, "Tissue Stratified" = tissue.stratified, "Tissue Specific" = res))
}

