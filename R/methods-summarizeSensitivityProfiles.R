#' Takes the sensitivity data from a PharmacoSet, and summarises them into a
#' drug vs cell line table
#'
#' This function creates a table with cell lines as rows and drugs as columns,
#' summarising the drug senstitivity data of a PharmacoSet into drug-cell line
#' pairs
#'
#' @examples
#' data(GDSCsmall)
#' GDSCauc <- summarizeSensitivityProfiles(GDSCsmall,
#'     sensitivity.measure='auc_published')
#'
#' @param object [PharmacoSet] The PharmacoSet from which to extract the data
#' @param sensitivity.measure `character` which sensitivity sensitivity.measure
#'   to use? Use the sensitivityMeasures function to find out what measures are
#'   available for each object.
#' @param cell.lines \code{character} The cell lines to be summarized.
#'   If any cell lines has no data, it will be filled with
#'   missing values
#' @param drugs \code{character} The drugs to be summarized.
#'   If any drugs has no data, it will be filled with
#'   missing values
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#'
#' @return [matrix] A matrix with cell lines going down the rows, drugs across
#'   the columns, with the selected sensitivity statistic for each pair.
#'
#' @importMethodsFrom CoreGx summarizeSensitivityProfiles
#' @export
setMethod("summarizeSensitivityProfiles", signature(object="PharmacoSet"),
        function(object, sensitivity.measure="auc_recomputed", cell.lines,
        drugs, summary.stat=c("mean", "median", "first", "last", "max", "min"),
        fill.missing=TRUE, verbose=TRUE) {
    if (is(treatmentResponse(object), 'LongTable'))
        .summarizeSensProfiles(object, sensitivity.measure,
            cell.lines, drugs, summary.stat, fill.missing)
    else
        .summarizeSensitivityProfilesPharmacoSet(object,
            sensitivity.measure, cell.lines, drugs, summary.stat,
            fill.missing, verbose)
})

#' Summarize the sensitivity profiles when the sensitivity slot is a LongTable
#'
#' @return [matrix] A matrix with cell lines going down the rows, drugs across
#'   the columns, with the selected sensitivity statistic for each pair.
#'
#' @import data.table
#' @keywords internal
.summarizeSensProfiles <- function(object,
        sensitivity.measure='auc_recomputed', cell.lines, drugs, summary.stat,
        fill.missing=TRUE) {

    # handle missing
    if (missing(cell.lines)) cell.lines <- sampleNames(object)
    if (missing(drugs)) drugs <- treatmentNames(object)
    if (missing(summary.stat)) summary.stat <- 'mean'

    # get LongTable object
    longTable <- treatmentResponse(object)

    # extract the sensitivty profiles
    sensProfiles <- assay(longTable, 'profiles', withDimnames=TRUE, key=FALSE)
    profileOpts <- setdiff(colnames(sensProfiles), idCols(longTable))

    # compute max concentration and add it to the profiles
    if (sensitivity.measure == 'max.conc') {
        dose <- copy(assay(longTable, 'dose', withDimnames=TRUE, key=FALSE))
        dose[, max.conc := max(.SD, na.rm=TRUE),
            .SDcols=grep('dose\\d+id', colnames(dose))]
        dose <- dose[, .SD, .SDcols=!grepl('dose\\d+id', colnames(dose))]
        sensProfiles <- dose[sensProfiles, on=idCols(longTable)]
    }

    # deal with drug combo methods
    if (sensitivity.measure == 'Synergy_score')
        drugs <- grep('///', drugs, value=TRUE)

    # ensure selected measure is an option
    if (!(sensitivity.measure %in% profileOpts))
        stop(.errorMsg('[PharmacoGx::summarizeSensivitiyProfiles,LongTable-method] ',
            'there is no measure ', sensitivity.measure, ' in this PharmacoSet.',
            ' Please select one of: ', .collapse(profileOpts)))

    # match summary function
    ## TODO:: extend this function to support passing in a custom summary function
    summary.function <- function(x) {
        if (all(is.na(x))) {
            return(NA_real_)
        }
        switch(summary.stat,
            "mean" = { mean(as.numeric(x), na.rm=TRUE) },
            "median" = { median(as.numeric(x), na.rm=TRUE) },
            "first" = { as.numeric(x)[[1]] },
            "last" = { as.numeric(x)[[length(x)]] },
            "max"= { max(as.numeric(x), na.rm=TRUE) },
            "min" = { min(as.numeric(x), na.rm=TRUE)}
            )
    }

    # do the summary
    profSummary <- sensProfiles[, summary.function(get(sensitivity.measure)),
        by=.(treatmentid, sampleid)]

    # NA pad the missing cells and drugs
    if (fill.missing) {
        allCombos <- data.table(expand.grid(drugs, cell.lines))
        colnames(allCombos) <- c("treatmentid", "sampleid")
        profSummary <- profSummary[allCombos, on=c("treatmentid", "sampleid")]
    }

    # reshape and convert to matrix
    setorderv(profSummary, c("sampleid", "treatmentid"))
    profSummary <- dcast(profSummary, treatmentid ~ sampleid, value.var='V1')
    summaryMatrix <- as.matrix(profSummary, rownames="treatmentid")
    return(summaryMatrix)

}



#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @importFrom reshape2 acast
#' @keywords internal
.summarizeSensitivityProfilesPharmacoSet <- function(object,
                                         sensitivity.measure="aac_recomputed",
                                         cell.lines,
                                         drugs,
                                         summary.stat=c("mean", "median", "first", "last", "max", "min"),
                                         fill.missing=TRUE, verbose=TRUE) {

	summary.stat <- match.arg(summary.stat)
  #sensitivity.measure <- match.arg(sensitivity.measure)
  if (!(sensitivity.measure %in% c(colnames(sensitivityProfiles(object)), "max.conc"))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", annotation(object)$name, paste(colnames(sensitivityProfiles(object)), collapse=", ")))
  }
  if (missing(cell.lines)) {
    cell.lines <- sampleNames(object)
  }
  if (missing(drugs)) {
    if (sensitivity.measure != "Synergy_score")
    {
      drugs <- treatmentNames(object)
    }else{
      drugs <- sensitivityInfo(object)[grep("///", sensitivityInfo(object)$treatmentid), "treatmentid"]
    }
  }

  pp <- sensitivityInfo(object)
  ppRows <- which(pp$sampleid %in% cell.lines & pp$treatmentid %in% drugs) ### NEEDED to deal with duplicated rownames!!!!!!!
  if(sensitivity.measure != "max.conc") {
    dd <- sensitivityProfiles(object)
  } else {

    if(!"max.conc" %in% colnames(sensitivityInfo(object))) {

      object <- updateMaxConc(object)

    }
    dd <- sensitivityInfo(object)

  }

  result <- matrix(NA_real_, nrow=length(drugs), ncol=length(cell.lines))
  rownames(result) <- drugs
  colnames(result) <- cell.lines

  # if(verbose){

  #   message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, annotation(object)$name))
  #   total <- length(drugs)*length(cell.lines)
  #   # create progress bar
  #   pb <- utils::txtProgressBar(min=0, max=total, style=3)
  #   i <- 1


  # }
  if(is.factor(dd[, sensitivity.measure]) | is.character(dd[, sensitivity.measure])){
    warning("Sensitivity measure is stored as a factor or character in the pSet. This is incorrect.\n
             Please correct this and/or file an issue. Fixing in the call of this function.")
    dd[, sensitivity.measure] <- as.numeric(as.character(dd[, sensitivity.measure]))
  }

  pp_dd <- cbind(pp[,c("sampleid", "treatmentid")], "sensitivity.measure"=dd[, sensitivity.measure])


  summary.function <- function(x) {
    if(all(is.na(x))){
      return(NA_real_)
    }
    switch(summary.stat,
        "mean" = { mean(as.numeric(x), na.rm=TRUE) },
        "median" = { median(as.numeric(x), na.rm=TRUE) },
        "first" = { as.numeric(x)[[1]] },
        "last" = { as.numeric(x)[[length(x)]] },
        "max"= { max(as.numeric(x), na.rm=TRUE) },
        "min" = { min(as.numeric(x), na.rm=TRUE)}
        )
  }

  pp_dd <- pp_dd[pp_dd[,"sampleid"] %in% cell.lines & pp_dd[,"treatmentid"]%in%drugs,]

  tt <- reshape2::acast(pp_dd, treatmentid ~ sampleid, fun.aggregate=summary.function, value.var="sensitivity.measure")
 # tt <- tt[drugs, cell.lines]



  result[rownames(tt), colnames(tt)] <- tt

	if (!fill.missing) {

    myRows <- apply(result, 1, function(x) !all(is.na(x)))
    myCols <- apply(result, 2, function(x) !all(is.na(x)))
    result <- result[myRows, myCols]
	}
  return(result)
}
