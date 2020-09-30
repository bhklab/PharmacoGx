# ==== PharmacoSet Class

#'`[`
#'
#' @examples
#' data(CCLEsmall)
#' CCLEsmall["WM1799", "Sorafenib"]
#'
#' @param x object
#' @param i Cell lines to keep in object
#' @param j Drugs to keep in object
#' @param ... further arguments
#' @param drop A boolean flag of whether to drop single dimensions or not
#'
#' @return Returns the subsetted object
#'
#' @export
setMethod(`[`, 'PharmacoSet', function(x, i, j, ..., drop = FALSE){
    eval(substitute(subset(x, i, j, i)))  # eval(substitute()) idiom allows
                                          # direct pass through of unevaluated
                                          # function arguments. Stops attempts
                                          # to evaluate arguments in current
                                          # scope.
})

# ==== LongTable Class

#' [ LongTable Method
#'
#' Single bracket subsetting for
#'
#' This function is endomorphic, it always returns a LongTable object.
#'
#' @param x [`LongTable`] The object to subset.
#' @param i [`character`], [`numeric`], [`logical`] or [`call`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all row id columns matching the vector. This parameter also
#'      supports valid R regex query strings which will match on the colnames
#'      of `x`. For convenience, * is converted to .* automatically. Colon
#'      can be to denote a specific part of the colnames string to query.
#'
#'  Numeric or Logical: these select based on the rowKey from the `rowData`
#'      method for the `LongTable`.
#'
#'  Call: Accepts valid query statements to the `data.table` i parameter as
#'      a call object. We have provided the function .() to conveniently
#'      convert raw R statements into a call for use in this function.
#'
#' @param j [`character`], [`numeric`], [`logical`] or [`call`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all drug id columns matching the vector. This parameter also
#'      supports regex queries. Colon can be to denote a specific part of the
#'      colnames string to query.
#'
#'  Numeric or Logical: these select based on the rowID from the `rowData`
#'      method for the `LongTable`.
#'
#'  Call: Accepts valid query statements to the `data.table` i parameter as
#'      a call object. We have provided the function .() to conveniently
#'      convert raw R statements into a call for use in this function.
#'
#' @export
setMethod('[', signature('LongTable'), function(x, i, j, assays, ..., drop=FALSE) {
    eval(substitute(subset(x, i, j, assays)))
})