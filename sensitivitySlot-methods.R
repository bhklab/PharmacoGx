#' sensitivitySlot Generic
#'
#' @examples
#' data(CCLEsmall)
#' sensitivitySlot(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#'
#' @return A \code{list} of the sensitivity slot contents
#'
#' @describeIn PharmacoSet Retrieve the contents of the sensitivity slot
#'
#' @importFrom methods callNextMethod
#' @importMethodsFrom CoreGx sensitivitySlot
#' @export
setMethod("sensitivitySlot", signature("PharmacoSet"), function(object) {
  callNextMethod(object)
})

##TODO:: Migrate this to CoreGx
#' sensitivitySlot<- Replacement Generic
#'
#' @examples
#' data(CCLEsmall)
#' sensitivitySlot(CCLEsmall) <- sensitivitySlot(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#' @param value A \code{list} of new sensitivity slot data for the pSet
#'
#' @return A copy of the \code{PharmacoSet} containing the updated sensitivty slot
#'
#' @describeIn PharmacoSet Set the raw dose and viability data for an pSet and return
#'   and updated copty
#'
#' @importFrom methods callNextMethod
#' @importMethodsFrom CoreGx sensitivitySlot<-
#' @export
setReplaceMethod("sensitivitySlot", signature("PharmacoSet", "list"),
                 function(object, value) {
                    callNextMethod(object=object, value=value)
                 })