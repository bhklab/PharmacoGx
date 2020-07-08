#' mDataNames Getter
#'
#' Returns the molecular data names for the `PharmacoSet` object
#'
#' @describeIn PharmacoSet Returns the names of molecular data types in a PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' mDataNames(CCLEsmall)
#'
#' @param object [`PharamcoSet`] object
#'
#' @return Vector of names of the molecular data types
#'
#' @importFrom CoreGx mDataNames
#' @importFrom methods callNextMethod
#'
#' @export
setMethod('mDataNames', signature='PharmacoSet', function(object) {
  callNextMethod(object=object)
})