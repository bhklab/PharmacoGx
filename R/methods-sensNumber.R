#' sensNumber Getter
#'
#' Get the sensitivity numbers for a `PharmacoSet` object
#'
#' @describeIn PharmacoSet Return the summary of available sensitivity
#'   experiments
#'
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall)
#'
#' @param object A \code{PharmacoSet}
#' @return A \code{data.frame} with the number of sensitivity experiments per
#'   drug and cell line
#'

#'
#' @importFrom CoreGx sensNumber
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(sensNumber, 'PharmacoSet', function(object){
  callNextMethod(object=object)
})