#' sensitivityRaw Getter
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityRaw(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#' @return A \code{array} containing the raw sensitivity data
#'
#' @describeIn PharmacoSet Retrive the raw dose and viability data from an pSet
#'
#' @importMethodsFrom CoreGx sensitivityRaw
#' @export
setMethod("sensitivityRaw", signature("PharmacoSet"), function(object) {
  callNextMethod(object)
})


#' sensitivityRaw<- Setter
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityRaw(CCLEsmall) <- sensitivityRaw(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#' @param value A \code{array} containing the raw dose and viability data for the
#'   pSet
#'
#' @return A copy of the \code{PharmacoSet} containing the updated sensitivty data
#'
#' @describeIn PharmacoSet Update the raw dose and viability data in an pSet object
#'
#' @importMethodsFrom CoreGx sensitivityRaw<-
#' @export
setReplaceMethod('sensitivityRaw', signature("PharmacoSet"), function(object, value) {
  callNextMethod(object=object, value=value)
})
