# ==== LongTable Class

#' [[ Method for LongTable Class
#'
#' Select an assay from within a LongTable object.
#'
#' @param x [`LongTable`] object to retrieve assays from
#' @param i [`character`] name or [`integer`] index of the desired assay.
#' @param keys [`logical`] Should the row and column keys also be returned?
#'    Defaults to FALSE to make applying aggregate functions more convenient.
#'
#' @describeIn LongTable
#'
#' @importFrom crayon cyan magenta
#' @export
setMethod('[[', signature('LongTable'), function(x, i, keys=FALSE) {
    if (length(i) > 1)
        stop(cyan$bold("Please only select one assay! To subset on multiple assays please
            see ?subset"))

    if (keys)
        return(assay(x, i))
    else
        return(assay(x, i)[, -c('rowKey', 'colKey')])
})