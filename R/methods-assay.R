#' Retrieve an assay `data.table` object from the `assays` slot of a `LongTable`
#'    object
#'
#' @param x [`LongTable`] The `LongTable` object to get the assay from.
#' @param i [`integer`] or [`character`] vector containing the index or name
#'   of the assay, respectively.
#' @param withDimnames [`logical`] Should the dimension names be returned
#'   joined to the assay. This retrieves both the row and column identifiers
#'   and returns them attached to the
#' @param metadata [`logical`] Should all of the metadata also be joined to
#'   the assay. This is useful when modifying assays as the resulting list
#'   has all the information needed to recreated the LongTable object.
#'
#' @importMethodsFrom SummarizedExperiment assay
#' @importFrom crayon magenta cyan
#' @export
setMethod('assay',
          signature(x='LongTable'),
          function(x, i, withDimnames=FALSE, metadata=FALSE) {

    # validate input
    if (length(i) > 1)
        stop(magenta$bold('Please specifying a single character assay name or
            integer index. See assayNames(x) for available assays.'))

    keepAssay <- if (is.character(i)) which(assayNames(x) == i) else i
    if (length(keepAssay) < 1)
        stop(magenta$bold('There is no assay named ',
                    i,
                    ' in this LongTable. Use assayNames(longTable) for a list of
                    valid assay names.'))

    # extract the specified assay
    assayData <- assays(x)[[keepAssay]]

    # optionally join to rowData and colData
    if (withDimnames && !metadata) {
        assayData <- .rowIDData(x)[assayData, on='rowKey'][, -'rowKey']
        assayData <- .colIDData(x)[assayData, on='colKey'][, -'colKey']
    } else if (withDimnames && metadata) {
        assayData <- rowData(x, key=TRUE)[assayData, on='rowKey'][, -'rowKey']
        assayData <- colData(x, key=TRUE)[assayData, on='colKey'][, -'colKey']
    }

    if (!withDimnames && metadata) warning(cyan$bold('Cannot use metadata=TRUE
        when withDimnames=FALSE. Ignoring the metadata argument.'))

    return(assayData)
})