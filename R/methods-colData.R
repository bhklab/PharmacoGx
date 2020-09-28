#' Retrieve the column metadata table from a LongTable object
#'
#' @param x A [`LongTable`] to retrieve column metadata from.
#' @param key [`logical`] She the colKey column also be returned? Defaults to
#'     FALSE.
#'
#' @return A [`data.table`] containing row identifiers and metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @import data.table
#' @export
setMethod('colData', signature(x='LongTable'), function(x, key=FALSE) {
    return(if (key) x@colData[, -'.colnames'] else
        x@colData[, -c('.colnames', 'colKey')])
})