#' Retrieve an assay `data.table` object from the `assays` slot of a `LongTable`
#'    object
#'
#' @param x [`LongTable`] The `LongTable` object to get the assay from.
#' @param i [`integer`] or [`character`] vector containing the index or name
#'   of the assay, respectively.
#' @param withDimnames [`logical`] Should the dimension names be returned
#'   joined to the assay. This retrieves both the row and column identifiers
#'   and returns them attached to the
#' @param
#'
#' @importMethodsFrom SummarizedExperiment assay
#' @export
setMethod('assay',
          signature(x='LongTable'),
          function(x, i, withDimnames=FALSE) {

    # validate input
    if (length(i) > 1)
        stop(paste0('Please specifying a single character assay name or
            integer index'))

    keepAssay <- if (is.character(i)) which(assayNames(x) == i) else i
    if (length(keepAssay) < 1)
        stop(paste0('There is no assay named ',
                    i,
                    ' in this LongTable. Use assayNames(longTable) for a list of
                    valid assay names.'))

    # extract the specified assay
    assayData <- assays(x)[[keepAssay]]

    # optionally join to rowData and colData
    if (withDimnames) {
        assayData <- .colIDData(x)[assayData, on='colKey'][, -'colKey']
        assayData <- .rowIDData(x)[assayData, on='rowKey'][, -'rowKey']
    }

    return(assayData)
})