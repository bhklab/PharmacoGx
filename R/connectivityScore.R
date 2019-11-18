########################
## Benjamin Haibe-Kains & Cat Olsen
## September 9, 2014
########################



#' Function computing connectivity scores between two signatures
#' 
#' A function for finding the connectivity between two signatures, using either
#' the GSEA method based on the KS statistic, or the gwc method based on a 
#' weighted spearman statistic. The GSEA analysis is implemented in the piano package. 
#' 
#' @references 
#'    F. Pozzi, T. Di Matteo, T. Aste, "Exponential smoothing weighted
#'    correlations", The European Physical Journal B, Vol. 85, No 6, 2012. DOI:
#'    10.1140/epjb/e2012-20697-x
#' @references
#'    Varemo, L., Nielsen, J. and Nookaew, I. (2013) Enriching the gene set
#'    analysis of genome-wide data by incorporating directionality of gene
#'    expression and combining statistical hypotheses and methods. Nucleic
#'    Acids Research. 41 (8), 4378-4391. doi: 10.1093/nar/gkt111
#'    
#' @examples
#' xValue <- c(1,5,23,4,8,9,2,19,11,12,13)
#' xSig <- c(0.01, 0.001, .97, 0.01,0.01,0.28,0.7,0.01,0.01,0.01,0.01)
#' yValue <- c(1,5,10,4,8,19,22,19,11,12,13)
#' ySig <- c(0.01, 0.001, .97,0.01, 0.01,0.78,0.9,0.01,0.01,0.01,0.01)
#' xx <- cbind(xValue, xSig)
#' yy <- cbind(yValue, ySig)
#' rownames(xx) <- rownames(yy) <- c('1','2','3','4','5','6','7','8','9','10','11')
#' data.cor <- connectivityScore(xx,yy,method="gwc", gwc.method="spearman", nperm=300)
#' 
#' @param x A \code{matrix} with the first gene signature. In the case of GSEA the vector of
#'   values per gene for GSEA in which we are looking for an enrichment. In the 
#'   case of gwc, this should be a matrix, with the per gene responses in the 
#'   first column, and the significance values in the second.
#' @param y A \code{matrix} with the second signature. In the case of GSEA, this is the
#'   vector of up and down regulated genes we are looking for in our signature,
#'   with the direction being determined from the sign. In the case of gwc, this
#'   should be a matrix of identical size to x, once again with the per gene
#'   responses in the first column, and their significance in the second.
#' @param method \code{character} string identifying which method to use, out of 'gsea' and 'gwc'
#' @param nperm \code{numeric}, how many permutations should be done to determine
#'   significance through permutation testing? The minimum is 100, default is
#'   1e4.
#' @param nthread \code{numeric}, how many cores to run parallel processing on.
#' @param gwc.method \code{character}, should gwc use a weighted spearman or pearson
#'   statistic?
#' @param ... Additional arguments passed down to gsea and gwc functions
#' @return \code{numeric} a numeric vector with the score and the p-value associated
#'   with it
#' @export
#' @importFrom piano runGSA
#' @importFrom piano loadGSC
#' @importFrom stats complete.cases



connectivityScore <- function(x, y, method=c("gsea", "fgsea", "gwc"), nperm=1e4, nthread=1, gwc.method=c("spearman", "pearson"), ...) {
  
  method <- match.arg(method)
  gwc.method <- match.arg(gwc.method)
  if (class(x) != "matrix") {
      x <- as.matrix(x)
  }
  if (class(y) != "matrix") {
      y <- as.matrix(y)
  }
  if ((ncol(x) != 2 || ncol(y) != 2) && method=="gwc") {
    stop ("x and y should have 2 columns: effect size and corresponding p-values")
  }
  if(method=="gsea"){
    method <- "fgsea"
    warning("Using fGSEA method to calculate GSEA")
  }
  if (method == "fgsea" && nrow(y) >= nrow(x)) {
    warning("GSEA method: query gene set (y) larger than signature (x)")
  }
  
  if (is.null(rownames(x)) || is.null(rownames(y)) || !length(intersect(rownames(x), rownames(y)))) {
    stop ("Row names of x and y are either missing or have no intersection")
  }
  if (nperm < 100){
    stop ("The minimum number of permutations for permutation testing is 100")
  }
  switch (method,
          "fgsea" = {
            ## remove missing values
            y <- y[!is.na(y[ ,1]), , drop=FALSE]
            x <- x[!is.na(x[ ,1]), , drop=FALSE]
            ## create gene set
            gset <- cbind("gene"=rownames(y), "set"=ifelse(as.numeric(y[ , 1]) >= 0, "UP", "DOWN")) 
            gset <- piano::loadGSC(gset)
            ## run enrichment analysis
            nes <- piano::runGSA(geneLevelStats=x[ , 1], geneSetStat="fgsea", gsc=gset, nPerm=nperm + (nperm %% nthread), ncpus=nthread, verbose=FALSE, adjMethod="none",...)
            ## merge p-values for negative and positive enrichment scores
            nes$pDistinctDir <- nes$pDistinctDirUp
            nes$pDistinctDir[is.na(nes$pDistinctDirUp), 1] <- nes$pDistinctDirDn[is.na(nes$pDistinctDirUp), 1]
            nes.up <- c(nes$statDistinctDir[which(names(nes$gsc) == "UP"), 1], nes$pDistinctDir[which(names(nes$gsc) == "UP"), 1])
            nes.down <- c(nes$statDistinctDir[which(names(nes$gsc) == "DOWN"), 1], nes$pDistinctDir[which(names(nes$gsc) == "DOWN"), 1])
            ## combine UP and DOWN
            if (length(nes.up) == 0){
              score = c("es" = -nes.down[1], "p" = nes.down[2])
            } else if (length(nes.down) == 0){
              score = c("es" = nes.up[1], "p" = nes.up[2])
            } else if (complete.cases(cbind(nes.up, nes.down)) && sign(nes.up[1]) != sign(nes.down[1])) {
              score <- c("es"=(nes.up[1] - nes.down[1]) / 2, "p"=combineTest(p=c(nes.up[2], nes.down[2]), method="fisher", na.rm=TRUE))
            } else {
              score <- c("score"=0, "p"=1)
            }
          },
          "gwc" = {
            ## intersection between x and y
            ii <- intersect(rownames(x), rownames(y))
            if(length(ii) < 10) {
              stop ("Less than 10 probes/genes in common between x and y")
            }
            score <- gwc(x1=x[ii, 1], p1=x[ii, 2], x2=y[ii, 1], p2=y[ii, 2], method.cor=gwc.method, nperm=nperm, ...)
            names(score) <- c("score", "p")
          }
  )
  return (score)
}
