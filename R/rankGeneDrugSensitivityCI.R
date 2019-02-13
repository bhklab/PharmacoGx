
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust

#################################################
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
#################################################

rankGeneDrugSensitivityCI <- function (data, drugpheno, verbose=FALSE) {


  if(is.null(dim(drugpheno))){

    drugpheno <- data.frame(drugpheno)

  } else if(class(drugpheno)!="data.frame"){
    drugpheno <- as.data.frame(drugpheno)

  }

  if (nrow(drugpheno) != nrow(data)) {
    stop("length of drugpheno should be equal to the number of rows of data!")
  }
  
  drugpheno <- drugpheno[[1]]
  res <- NULL
  # browser()
  res <- apply(data, 2, function(gene.exp) {
    return(justFastCI(gene.exp, drugpheno, noise.ties=TRUE))
  })

}

## End
