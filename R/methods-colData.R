#' Retrieve the column metadata table from a LongTable object
#'
#' @param x [`LongTable`]
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @import data.table
#' @export
setMethod('colData', signature(x='LongTable'), function(x) {
    return(x@colData[, -'.colnames'])
})