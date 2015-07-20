########################
## Petr Smirnov
## All rights Reserved
## September 3. 2014
########################

#	Get Drug Perturbation Signatures from a PharmacoSet
###############################################################################
## Drug perturbation analysis
## create profiles before vs after drug for each drug 
###############################################################################



#' Creates a signature representing gene expression (or other molecular profile)
#' change induced by administrating a drug, for use in drug effect analysis.
#' 
#' Given a Pharmacoset of the perturbation experiment type, and a list of drugs,
#' the function will compute a signature for the effect of drug concentration on
#' the molecular profile of a cell. The algorithm uses a regression model which
#' corrects for experimental batch effects, cell specific differences, and
#' duration of experiment to isolate the effect of the concentration of the drug
#' applied. The function returns the estimated coefficient for concentration,
#' the t-stat, the p-value and the false discovery rate associated with that
#' coefficient, in a 3 dimensional array, with genes in the first direction,
#' drugs in the second, and the selected return values in the third.
#' 
#' @examples
#' data(CMAPsmall)
#' drug.perturbation <- drugPertubrationSig(CMAPsmall)
#' head(drug.perturbation)
#' 
#' @param pSet \code{PharmacoSet} of the perturbation experiment type
#' @param molecularData \code{character} string, which one of the molecular data
#'   types to use in the analysis, out of DNA, RNA, SNP, CNV
#' @param drugs \code{character} vector of drug names for which to compute the 
#'   signatures. Should match the names used in the PharmacoSet.
#' @param nbcore \code{numeric}, if multiple cores are available, how many cores
#'   should the computation be parallelized over?
#' @param returnValues \code{character} vector, identifying which of estimate,
#'   t-stat, p-value and fdr should the function return for each gene drug pair?
#' @param verbose \code{bool} Should the function print diagnostic messages?
#'   Defaults to FALSE
#' @return A 3D \code{array} with genes in the first dimension, drugs in the 
#'   second, and return values in the third.
#' @export
#' @import parallel

drugPertubrationSig <- function(pSet, molecularData=c("rna", "dna", "snp", "cnv"), drugs, nbcore=1, returnValues=c("estimate","tstat", "pvalue", "fdr"), verbose=FALSE){
  
  molecularData <- match.arg(molecularData)
  
  availcore <- parallel::detectCores()
  if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
  options("mc.cores"=nbcore)
  
  if(missing(drugs)){
    drugn <- drugNames(pSet)
  } else {
    drugn <- drugs
  }
  
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
  
	dix <- is.element(drugn, pData(eset)[ , "drugid"])
  if (verbose && !all(dix)) {
    warning (sprintf("%i/%i drugs can be found", sum(dix), length(drugn)))
  }
  drugn <- drugn[dix]
  
  splitix <- parallel::splitIndices(nx=length(drugn), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, drugn, exprs, sampleinfo) {
    res <- NULL
    for(i in drugn[x]) {
      ## using a linear model (x ~ concentration + cell + batch + duration)
      rr <- rankGeneDrugPerturbation(data=exprs, drug=i, drug.id=as.character(sampleinfo[ , "drugid"]), drug.concentration=as.numeric(sampleinfo[ , "concentration"]), type=as.character(sampleinfo[ , "cellid"]), xp=as.character(sampleinfo[ , "xptype"]), batch=as.character(sampleinfo[ , "batchid"]), duration=as.character(sampleinfo[ , "duration"]) ,single.type=FALSE, nthread=1, verbose=FALSE)$all[ , returnValues, drop=FALSE]
      res <- c(res, list(rr))
    }
    names(res) <- drugn[x]
    return(res)
  }, drugn=drugn, exprs=t(exprs(eset)), sampleinfo=pData(eset))
  tt <- do.call(c, mcres)
  tt <- tt[!sapply(tt, is.null)]
  drug.perturbation <- array(NA, dim=c(nrow(geneInfo(pSet)), length(tt), ncol(tt[[1]])), dimnames=list(rownames(geneInfo(pSet)), names(tt), colnames(tt[[1]])))
  for(j in 1:ncol(tt[[1]])) {
    ttt <- sapply(tt, function(x, j, k) {
      xx <- array(NA, dim=length(k), dimnames=list(k))
      xx[rownames(x)] <- x[ , j]
      return (xx)
    }, j=j, k=rownames(geneInfo(pSet)))
    drug.perturbation[rownames(ttt), colnames(ttt), j] <- ttt
  }
  
  drug.perturbation <- PharmacoGxSignatures(drug.perturbation, PSetName = pSetName(pSet), Call = as.character(match.call()), SigType='Perturbation')
  
  return(drug.perturbation)
}
