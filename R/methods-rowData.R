#' Retrieve the row metadata table from a LongTable object
#'
#' @param x A [`LongTable`] object to retrieve the row metadata from.
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @import data.table
#' @export
setMethod('rowData', signature(x='LongTable'), function(x) {
    return(x@rowData[, -'.rownames'])
})