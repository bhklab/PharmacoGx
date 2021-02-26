#' sensitivityMeasures Getter
#'
#' Get the types of sensitivity measurements from a object object
#'
#' @examples
#' data(CCLEsmall)
#' sensMeas <- sensitivityMeasures(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] object to retrieve available sensitivity measurement types from
#'
#' @return A [`character`] vector of all the available sensitivity measures
#'
#' @importMethodsFrom CoreGx sensitivityMeasures
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#' @export
setMethod('sensitivityMeasures', "PharmacoSet", function(object){
  colnames(sensitivityProfiles(object))
})

#' sensitivityMeasures Setter
#'
#' Get the types of sensitivity measurements available in a PharmacoSet object
#'
#' @describeIn PharmacoSet returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#'
#' @examples
#' data(CCLEsmall)
#' sensMeas <- sensitivityMeasures(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] object to retrieve available sensitivity measurement types from
#'
#' @return A [`character`] vector of all the available sensitivity measures
#'
#' @importMethodsFrom CoreGx sensitivityMeasures<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod('sensitivityMeasures',
                 signature=signature(object='PharmacoSet', value='character'),
                 function(object, value) {
  colnames(sensitivityProfiles(object)) <- value
})