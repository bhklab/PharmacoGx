# ==== LongTable Class

#' Getter method for the metadata slot of a `LongTable` object
#'
#' @param x The [`LongTable`] object from which to retrieve the metadata list.
#'
#' @return [`list`] The contents of the `metadata` slot of the `LongTable`
#'   object.
#'
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata', signature(x='LongTable'), function(x) {
    return(x@metadata)
})

#' Setter method for the metadata slot of a `LongTable` object
#'
#' @param x [`LongTable`] The LongTable to update
#' @param value [`list`] A list of new metadata associated with a `LongTable`
#'   object.
#'
#' @return [`LongTable`] A copy of the `LongTable` object with the `value` in
#'   the metadata slot.
#'
#' @importMethodsFrom S4Vectors metadata<-
#' @importFrom crayon cyan magenta
#' @export
setReplaceMethod('metadata', signature(x='LongTable'), function(x, value) {
    if (!is(value, 'list'))
        stop(magenta$bold('The `metadata` slot must be a list!'))
    x@metadata <- value
    return(x)
})