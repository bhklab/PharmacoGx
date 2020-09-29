# ==== LongTable Class

#' Get the column names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to get the column names from
#'
#' @return [`character`] Vector of column names.
#'
#' @export
setMethod('colnames', signature(x='LongTable'), function(x) {
    return(x@colData$.colnames)
})

#' Get the row names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to get the row names from
#'
#' @return [`character`] Vector of row names.
#'
#' @export
setMethod('rownames', signature(x='LongTable'), function(x) {
    return(x@rowData$.rownames)
})

#' Getter for the dimnames of a `LongTable` object
#'
#' @param x The [`LongTable`] object to retrieve the dimnames for
#'
#' @return [`list`] List with two character vectors, one for row and one for
#'     column names.
#'
#' @importMethodsFrom Biobase dimnames
#' @export
setMethod('dimnames', signature(x='LongTable'), function(x) {
    return(list(rownames(x), colnames(x)))
})