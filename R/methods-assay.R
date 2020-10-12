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
#' @param key [`logical`] Should the key columns also be returned?
#'
#' @importMethodsFrom SummarizedExperiment assay
#' @importFrom crayon magenta cyan
#' @export
##TODO:: Add key argument with default to FALSE to remove rowKey and colKey
setMethod('assay',
          signature(x='LongTable'),
          function(x, i, withDimnames=FALSE, metadata=FALSE, key=TRUE) {

    # validate input
    if (length(i) > 1)
        stop(.errorMsg('\nPlease specifying a single character assay name or',
            'integer index. See assayNames(x) for available assays.'))

    keepAssay <- if (is.character(i)) which(assayNames(x) == i) else i
    if (length(keepAssay) < 1)
        stop(.errorMsg('\nThere is no assay ', i,
            ' in this LongTable. Use assayNames(longTable) for a list',
            'of valid assay names.'))

    # extract the specified assay
    assayData <- x@assays[[keepAssay]]

    # optionally join to rowData and colData
    if (withDimnames && !metadata) {
        assayData <- rowIDs(x, data=TRUE, key=TRUE)[assayData, on='rowKey']
        assayData <- colIDs(x, data=TRUE, key=TRUE)[assayData, on='colKey']
    } else if (withDimnames && metadata) {
        assayData <- rowData(x, key=TRUE)[assayData, on='rowKey']
        assayData <- colData(x, key=TRUE)[assayData, on='colKey']
    }

    if (!key) assayData <- assayData[, -c('rowKey', 'colKey')]

    if (!withDimnames && metadata)
    warning(.warnMsg('\nCannot use metadata=TRUE when withDimnames=FALSE.',
        'Ignoring the metadata argument.'))

    assayData

})


#' Add or replace an assay in a LongTable by name
#'
#' @param x [`LongTable`]
#' @param i [`character`]
#' @param value [`data.frame`] or [`data.table`]
#'
#' @return [`LongTable`] With updated assays slot.
#'
#' @importMethodsFrom SummarizedExperiment assay<-
#' @export
setReplaceMethod('assay',
                 signature(x='LongTable', i='character'),
                 function(x, i, value) {

    if (!is.data.frame(value)) stop(.errorMsg('Only a data.frame or data.table
        can be assiged to the assay slot.'))

    if (length(i) > 1) stop(.errorMsg("Only a single assay name can be assiged
        with assay(x, i) <- value."))

    whichAssay <- which(i %in% assayNames(x))

    assayData <- assays(x, withDimnames=TRUE, metadata=TRUE)

    if (!is.data.table(value)) setDT(value)

    if (length(whichAssay) > 0) {
        assayData[[i]] <- value
    } else {
        assayData <- c(assayData, eval(str2lang(paste0('list(', i, '=value)'))))
    }

    ## TODO:: Do I need to call the build function?

    x@assays <- assayData
    x
})