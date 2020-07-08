#' featureInfo getter
#'
#' Get the molecular profile data for the specified molecular data type
#'
#' @examples
#' data(CCLEsmall)
#' featInf <- featureInfo(CCLEsmall, "rna")
#'
#' @param object The [`PharmacoSet`] to retrieve feature annotations from
#' @param mDataType [`character`] The type of molecular data
#' @return A [`data.frame`] with the feature annotations for the specified `mDataType`
#'
#' @describeIn PharmacoSet Return the feature info for the given molecular datatype
#'
#' @importFrom methods callNextMethod
#' @importMethodsFrom CoreGx featureInfo
#' @export
setMethod(featureInfo,
          signature="PharmacoSet",
          function(object, mDataType){
  callNextMethod(object=object, mDataType=mDataType)
})

#' featureInfo<- setter
#'
#' Update the molecular profile data for the specified datatype in the specified pSet object
#'
#' @examples
#' data(CCLEsmall)
#' featureInfo(CCLEsmall, "rna") <- featureInfo(CCLEsmall, "rna")
#'
#' @param object The [`PharmacoSet`] to replace gene annotations in
#' @param mDataType [`character`] The type of molecular data to be updated
#' @param value A [`DataFrame`] with the new feature annotations
#'
#' @return Updated \code{PharmacoSet}
#'
#' @importMethodsFrom CoreGx featureInfo<-
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Replace the gene info for the molecular data
#'
#' @export
setReplaceMethod("featureInfo",
                 signature=signature(object="PharmacoSet", mDataType ="character", value="DataFrame"),
                 function(object, mDataType, value) {
  callNextMethod(object=object, mDataType=mDataType, value=value)
})
