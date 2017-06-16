#' Takes the sensitivity data from a PharmacoSet, and summarises them into a
#' drug vs cell line table
#' 
#' This function creates a table with cell lines as rows and drugs as columns,
#' summarising the drug senstitivity data of a PharmacoSet into drug-cell line
#' pairs
#' 
#' @examples 
#' data(GDSCsmall)
#' GDSCauc <- summarizeSensitivityProfiles(GDSCsmall, sensitivity.measure='auc_published')
#'
#' @param pSet [PharmacoSet] The PharmacoSet from which to extract the data
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? Use the 
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param cell.lines \code{character} The cell lines to be summarized. 
#'    If any cell lines has no data, it will be filled with
#'   missing values
#' @param drugs \code{character} The drugs to be summarized.
#'   If any drugs has no data, it will be filled with
#'   missing values
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#' @return [matrix] A matrix with cell lines going down the rows, drugs across
#'   the columns, with the selected sensitivity statistic for each pair.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @importFrom reshape2 acast
#' @export


summarizeSensitivityProfiles <- function(pSet, sensitivity.measure="auc_recomputed", cell.lines, drugs, summary.stat=c("mean", "median", "first", "last", "max", "min"), fill.missing=TRUE, verbose=TRUE){
	summary.stat <- match.arg(summary.stat)
  #sensitivity.measure <- match.arg(sensitivity.measure)
  if (!(sensitivity.measure %in% c(colnames(sensitivityProfiles(pSet)),"max.conc"))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", pSet@annotation$name, paste(colnames(sensitivityProfiles(pSet)), collapse=", ")))
  }
  if (missing(cell.lines)) {
    cell.lines <- cellNames(pSet)
  }
  if (missing(drugs)) {
    if (sensitivity.measure != "Synergy_score")
    {
      drugs <- drugNames(pSet)
    }else{
      drugs <- sensitivityInfo(pSet)[grep("///", sensitivityInfo(pSet)$drugid), "drugid"]
    }
  }
  
  pp <- sensitivityInfo(pSet)
  ppRows <- which(pp$cellid %in% cell.lines & pp$drugid %in% drugs) ### NEEDED to deal with duplicated rownames!!!!!!!
  if(sensitivity.measure != "max.conc") {
    dd <- sensitivityProfiles(pSet)
  } else {

    if(!"max.conc"%in% colnames(sensitivityInfo(pSet))){

      pSet <- updateMaxConc(pSet)

    }
    dd <- sensitivityInfo(pSet)

  }

  result <- matrix(NA_real_, nrow=length(drugs), ncol=length(cell.lines))
  rownames(result) <- drugs
  colnames(result) <- cell.lines

  # if(verbose){

  #   message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, pSet@annotation$name))
  #   total <- length(drugs)*length(cell.lines)
  #   # create progress bar 
  #   pb <- utils::txtProgressBar(min=0, max=total, style=3)
  #   i <- 1


  # }

  pp_dd <- cbind(pp[,c("cellid", "drugid")], "sensitivity.measure"=dd[, sensitivity.measure])


  summary.function <- function(x) {
    if(all(is.na(x))){
      return(NA_real_)
    }
    switch(summary.stat, 
        "mean" = {
          return(mean(as.numeric(x), na.rm=TRUE))
        },
        "median" = {
          return(median(as.numeric(x), na.rm=TRUE))
        }, 
        "first" = {
          return(as.numeric(x)[[1]])
        },
        "last" = {
          return(as.numeric(x)[[length(x)]])
        },
        "max"= {
          return(max(as.numeric(x), na.rm=TRUE))
        },
        "min" = { 
          return(min(as.numeric(x), na.rm=TRUE))
        })

  }
  
  pp_dd <- pp_dd[pp_dd[,"cellid"]%in%cell.lines & pp_dd[,"drugid"]%in%drugs,]

  tt <- reshape2::acast(pp_dd, drugid~cellid, fun.aggregate=summary.function, value.var="sensitivity.measure")
 # tt <- tt[drugs, cell.lines]
  
  

  result[rownames(tt), colnames(tt)] <- tt

	if (!fill.missing) {
	  
    myRows <- apply(result, 1, function(x) !all(is.na(x)))
    myCols <- apply(result, 2, function(x) !all(is.na(x)))
    result <- result[myRows, myCols]
	}
  return(result)
}
