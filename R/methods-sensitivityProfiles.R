#' sensitivityProfiles getter
#'
#' Get the sensitivityProfiles data.frame from a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensProf <- sensitivityProfiles(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to retrieve sensitivity experiment data from
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @importFrom CoreGx sensitivityProfiles
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the phenotypic data for the drug dose sensitivity
#' @export
setMethod(sensitivityProfiles, "PharmacoSet", function(object) {
  callNextMethod(object=object)
})


#' sensitivityProfiles<- setter
#'
#' Set the sensitivityProfiles in a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityProfiles(CCLEsmall) <- sensitivityProfiles(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{data.frame} with the new sensitivity profiles. If a matrix object is passed in, converted to data.frame before assignment
#'
#' @return Updated \code{PharmacoSet}
#'
#' @importFrom CoreGx sensitivityProfiles<-
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet", value="data.frame"), function(object, value){
  callNextMethod(object=object, value=value)
})
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet", value="matrix"), function(object, value){
  callNextMethod(object=object, value=value)
})