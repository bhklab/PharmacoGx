#' dateCreated Getter
#'
#' Get the data that a `PharmacoSet` object was updated
#' @describeIn PharmacoSet Return the date the PharmacoSet was created
#'
#' @examples
#' data(CCLEsmall)
#' dateCreated(CCLEsmall)
#'
#' @param object A [`PharmacoSet`]
#'
#' @return [`character`] The date the `PharmacoSet` was created
#'
#' @importMethodsFrom CoreGx dateCreated
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(dateCreated, 'PharmacoSet', function(object) {
  callNextMethod(object=object)
})


#' sensNumber<- Generic
#'
#' A generic for the sensNumber method
#'
#'
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall) <- sensNumber(CCLEsmall)
#'
#' @param object A \code{PharmacoSet}
#' @param value A new \code{data.frame} with the number of sensitivity experiments per drug and cell line
#'
#' @return The updated \code{PharmacoSet}
#'
#' @describeIn PharmacoSet Update the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod('sensNumber', signature = signature(object='PharmacoSet',value='matrix'), function(object, value){
  callNextMethod(object=object, value=value)
})