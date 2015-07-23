## TODO::: NEED TO DEAL WITH DUPLICATES 
## TODO:: SOMEONE WITH MORE KNOWLEDGE OF
## THESE FUNCTIONS SHOULD WRITE BETTER DOCUMENTATION

########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:	
##      - data: gene expression data matrix
##			- drugpheno: sensititivity values fo thr drug of interest
##			- type: cell or tissue type for each experiment
##			- duration: experiment duration in hours
##      - batch: experiment batches
##			- single.type: Should the statitsics be computed for each cell/tissue type separately?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs:
## list of datafraes with the statistics for each gene, for each type
##
## Notes:	duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
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
#' 
#' @param pSet \code{PharmacoSet} of the perturbation experiment type
#' @param molecularData \code{character} string, which one of the molecular data
#'   types to use in the analysis, out of DNA, RNA, SNP, CNV
#' @param drugs \code{character} vector of drug names for which to compute the 
#'   signatures. Should match the names used in the PharmacoSet.
#' @param nbcore \code{numeric}, if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues \code{character} vector, detailing which of estimate,
#'   t-stat, p-value and fdr should the function return for each gene drug pair?
#' @param sensitivity.measure \code{character} string, idetifying which measure
#'   of the drug dose sensitivity should the function use for its computations?
#'   The current choices are 'ic50_published', 'auc_published',
#'   'ic50_recomputed', 'auc_recomputed'.
#' @param duplicates \code{character} string, identifying which summary
#'   statistic should be used to summarize duplicates for cell line sensitivity
#'   measurements? Currently implemented is only the median.
#' @param verbose \code{bool} Should the function print diagnostic messages?
#' @return A 3D \code{array} with genes in the first dimension, drugs in the 
#'   second, and return values in the third.
#' @export
#' @import parallel

drugSensitivitySig <- function(pSet, molecularData=c("rna", "dna", "snp", "cnv"), drugs, sensitivity.measure=c("ic50_published", "auc_published", "ic50_recomputed", "auc_recomputed"), duplicates="median", returnValues=c("estimate","tstat", "pvalue", "fdr"), nbcore=1, verbose=FALSE) {
	
	### This function needs to: Get a table of AUC values per cell line / drug
	### Be able to recompute those values on the fly from raw data if needed to change concentration
	### Be able to choose different summary methods on fly if needed (need to add annotation to table to tell what summary method previously used)
	### Be able to extract genomic data 
	### Run rankGeneDrugSens in parallel at the drug level
	### Return matrix as we had before
	
  
  molecularData <- match.arg(molecularData)
  sensitivity.measure <- match.arg(sensitivity.measure)
  
  if(missing(drugs)){
    drugn <- drugNames(pSet)
  } else {
    drugn <- drugs
  }
  
	drugpheno.all <- summarizeSensitivityPhenotype(pSet, sensitivity.measure=sensitivity.measure)
	dix <- is.element(drugn, colnames(drugpheno.all))
  if (verbose && !all(dix)) {
    warning (sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
  }
  drugn <- drugn[dix]

	switch (molecularData, 
    "rna" = {
      eset <- pSet@molecularData$rna
      },
    "dna" = {
      stop ("Drug sensitivity signature for DNA is not implemented yet")
      eset <- pSet@molecularData$dna
      },
    "snp" = {
      stop ("Drug sensitivity signature for SNP is not implemented yet")
      eset <- pSet@molecularData$snp
      },
    "cnv" = {
      stop ("Drug sensitivity signature for CNV is not implemented yet")
      eset <- pSet@molecularData$cnv
    }
  )
  
	if (any(duplicated(pData(eset)[,"cellid"]))){
		
		dupl_cells = unique(pData(eset)[which(duplicated(pData(eset)[,"cellid"])),"cellid"])
		for (cell in dupl_cells){
			myx <- which(pData(eset)[,"cellid"]==cell)
			med_exprs <- median(exprs(eset)[,myx])
			exprs(eset)[,myx[1]] <- med_exprs 
		}
		exprs(eset) <- exprs(eset)[,!duplicated(pData(eset)[,"cellid"])]
		pData(eset) <- pData(eset)[!duplicated(pData(eset)[,"cellid"]),]
	}
	pSet@molecularData$rna <- eset
	
	drugpheno.all <- drugpheno.all[rnaInfo(pSet)[,"cellid"],, drop=FALSE]
	
	type <- as.factor(cellInfo(pSet)[rnaInfo(pSet)[,"cellid"], "tissue.type"]) 
    batch <- rnaInfo(pSet)[, "batch"]
    batch[!is.na(batch) & batch == "NA"] <- NA
    batch <- as.factor(batch)
	names(batch) <- rnaInfo(pSet)[ , "cellid"]
	batch <- batch[rownames(drugpheno.all)]
	# duration <- sensitivityInfo(pSet)[,"duration_h"]
    ## compute drug sensitivity signatures
    availcore <- parallel::detectCores()
    if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
    options("mc.cores"=nbcore)
    message("Compute sensitivity signatures for all drugs")
    splitix <- parallel::splitIndices(nx=length(drugn), ncl=nbcore)
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <-  parallel::mclapply(splitix, function(x, drugn, expr, drugpheno, type, batch) {
      res <- NULL
      for(i in drugn[x]) {
        ## using a linear model (x ~ concentration + cell + batch)
        rr <- rankGeneDrugSensitivity(data=expr, drugpheno=drugpheno[ , i], type=type, batch=batch, single.type=FALSE, nthread=1)
        res <- c(res, list(rr$all))
      }
      names(res) <- drugn[x]
      return(res)
    }, drugn=drugn, expr=t(rnaData(pSet)), drugpheno=drugpheno.all, type=type, batch=batch)
    res <- do.call(c, mcres)
    res <- res[!sapply(res, is.null)]
    drug.sensitivity <- array(NA, dim=c(nrow(geneInfo(pSet)), length(res), ncol(res[[1]])), dimnames=list(rownames(geneInfo(pSet)), names(res), colnames(res[[1]])))
    for(j in 1:ncol(res[[1]])) {
      ttt <- sapply(res, function(x, j, k) {
        xx <- array(NA, dim=length(k), dimnames=list(k))
        xx[rownames(x)] <- x[ , j]
        return (xx)
      }, j=j, k=rownames(geneInfo(pSet)))
      drug.sensitivity[rownames(ttt), colnames(ttt), j] <- ttt
    }

drug.perturbation <- PharmacoSig(drug.perturbation, PSetName = pSetName(pSet), Call =match.call(), SigType='Sensitivity')

	return(drug.sensitivity)
}
