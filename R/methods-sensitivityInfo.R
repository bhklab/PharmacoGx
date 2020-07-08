#' sensitivityInfo getter
#'
#' Get the senstivity information DataFrame from a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensInf <- sensitivityInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to retrieve sensitivity experiment annotations from
#'
#' @return a [`DataFrame`] with the experiment info
#'
#' @importMethodsFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the drug dose sensitivity experiment info
#'
#' @export
setMethod(sensitivityInfo, signature="PharmacoSet", function(object){
  callNextMethod(object=object)
})


#' sensitivityInfo<- setter
#'
#' Set the sensitivityInfo DataFrame in a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{DataFrame} with the new sensitivity annotations
#'
#' @return Updated \code{PharmacoSet}
#'
#' @importFrom CoreGx sensitivityInfo<-
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Update the sensitivity experiment info
#' @export
setReplaceMethod("sensitivityInfo",
                 signature = signature(object="PharmacoSet",value="data.frame"),
                 function(object, value){
  callNextMethod(object=object, value=value)
})