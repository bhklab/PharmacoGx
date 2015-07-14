########################
## Petr Smirnov
## All rights Reserved
## June 12, 2015
########################
### MAYBE WE NEED A MORE USEFULL EXAMPLE HERE

##' Takes the gene expression RNA data from a PharmacoSet, and summarises them
##' into one entry per drug
##' 
##' Given a PharmacoSet with RNA expression data, this function will summarize
##' the data into one profile per cell line, using the chosed summaryStat. Note
##' that this does not really make sense with perturbation type data, and will
##' combine experiments and controls when doing the summary if run on a
##' perturbation dataset.
##' 
##' @examples
##' data(CGPsmall)
##' CGPsmall <- summarizeGeneExpression(CGPsmall, cells=cellNames(CGPsmall),
##'                     summaryStat = 'median', fillMissing = TRUE, verbose=TRUE)
##' CGPsmall
##'
##' @param pSet \code{PharmacoSet} The PharmacoSet to summarize
##' @param cells \code{character} The cell lines to be summarized in the 
##'   resulting PharmacoSet. If any cells have no data, they will be filled with
##'   missing values
##' @param summaryStat \code{character} which summary method to use if there are repeated
##'   cells?
##' @param fillMissing \code{boolean} should the missing cell lines not in the
##'   expression object be filled in with missing values?
##' @param verbose \code{bool} Should the function report progress?
##' @return \code{matrix} An updated PharmacoSet with the RNA expression summarized
##'   per cell line.
##' @export



summarizeGeneExpression <- function(pSet, cells, summaryStat="median", fillMissing=TRUE, verbose=TRUE){
  
  sampleinfo <- rnaInfo(pSet)[vector("numeric"),]
  sampleData <- rnaData(pSet)[,vector("numeric")]
  if (!is.null(cells)){
    if (verbose){
      cat("Summarising Gene Expression for:\t", pSet@annotation$name, "\n")
      total <- length(cells)
      # create progress bar 
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      i = 1
    }
    for (x in cells) {
      
      myx <- which(rnaInfo(pSet)[,"cellid"]==x)
      
      if (length(myx)==1){
        cellExprs <- rnaData(pSet)[,myx, drop=FALSE]
        cellAnnot <- rnaInfo(pSet)[myx[1],, drop=FALSE]
        sampleinfo <- rbind(sampleinfo, cellAnnot)
        sampleData <- cbind(sampleData, cellExprs)
      }
      if (length(myx)==0){
        if (fillMissing){
          cellExprs <- rep(NA, times=nrow(rnaData(pSet)))
          cellAnnot <- matrix(NA, nrow=1, ncol=ncol(rnaInfo(pSet)), dimnames=list(c(x), colnames=colnames(rnaInfo(pSet))))
          cellAnnot[,"cellid"] <- x
          sampleinfo <- rbind(sampleinfo, cellAnnot)
          sampleData <- cbind(sampleData, cellExprs)
        }
      } else if (length(myx)>1) {
        switch(summaryStat, "mean"={
          
          cellExprs <- apply(rnaData(pSet)[,myx, drop=FALSE], 1, mean)
          cellAnnot <- rnaInfo(pSet)[myx[1],,drop=FALSE]
          sampleinfo <- rbind(sampleinfo, cellAnnot)
          sampleData <- cbind(sampleData, cellExprs)
        },
        "median"={
          
          cellExprs <- apply(rnaData(pSet)[,myx, drop=FALSE], 1, median)
          cellAnnot <- rnaInfo(pSet)[myx[1],,drop=FALSE]
          sampleinfo <- rbind(sampleinfo, cellAnnot)
          sampleData <- cbind(sampleData, cellExprs)
        }, 
        "first"={
          
          cellExprs <- rnaData(pSet)[,myx[1],drop=FALSE]
          cellAnnot <- rnaInfo(pSet)[myx[1],,drop=FALSE]
          sampleinfo <- rbind(sampleinfo, cellAnnot)
          sampleData <- cbind(sampleData, cellExprs)
        },
        "last" = {
          
          cellExprs <- rnaData(pSet)[,myx[length(myx)], drop=FALSE]
          cellAnnot <- rnaInfo(pSet)[myx[length(myx)],,drop=FALSE]
          sampleinfo <- rbind(sampleinfo, cellAnnot)
          sampleData <- cbind(sampleData, cellExprs)
        })
        
      }
      if (verbose){
        setTxtProgressBar(pb, i)
        i = i+1
      }
    }
  }
  colnames(sampleData) <- cells
  rownames(sampleinfo) <- cells
  rnaData(pSet) <- sampleData
  rnaInfo(pSet) <- sampleinfo
  if (verbose) {close(pb)}
  return(pSet)
}

