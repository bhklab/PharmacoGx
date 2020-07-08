#' cellNames Getter
#'
#' Get the names of all cell-lines available in a `PharmacoSet` object
#'
#' @describeIn PharmacoSet Return the cell names used in the dataset
#'
#' @examples
#' data(CCLEsmall)
#' cellNames(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to return cell names from
#' @return A vector of the cell names used in the PharmacoSet
#'
#' @importMethodsFrom CoreGx cellNames
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(cellNames, signature="PharmacoSet", function(object) {
  callNextMethod(object)
})


#' cellNames<- Setter
#'
#' Update the names of cell lines available in a `PharmacoSet` object
#'
#' @examples
#' data(CCLEsmall)
#' cellNames(CCLEsmall) <- cellNames(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] object to update
#' @param value A [`character`] vector of the new cell names
#'
#' @return Updated [`PharmacoSet`]
#'
#' @importMethodsFrom CoreGx cellNames<-
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Update the cell names used in the dataset
#' @export
setReplaceMethod("cellNames",
                 signature=signature(object="PharmacoSet",value="character"),
                 function(object, value) {
    callNextMethod(object, value=value)
})