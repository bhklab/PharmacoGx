#' fNames
#'
#' Return the feature names for the specified molecular data type
#'
#' @describeIn PharmacoSet Return the feature names used in the dataset
#'
#' @examples
#' data(CCLEsmall)
#' fNames(CCLEsmall, "rna")
#'
#' @param object The [`PharmacoSet`]
#' @param mDataType [`character`] The molecular data type to return feature names for
#'
#' @return A [`character`] vector of the feature names
#'
#' @importFrom CoreGx fNames
#' @importFrom methods callNextMethod
#'
#' @export
setMethod('fNames',
          signature=signature(object='PharmacoSet', mDataType='character'),
          function(object, mDataType){
    callNextMethod(object=object, mDataType=mDataType)
})


#' fNames<-
#'
#' Setter for the feature names of a [`SummarizedExperiment`] in the
#'   molecularProfiles slot
#'
#' @examples
#' data(CCLEsmall)
#' fNames(CCLEsmall, 'rna') <- fNames(CCLEsmall, 'rna')
#'
#' @param object The [`PharmacoSet`] object to update
#' @param mDataType [`character`] The molecular data type to update
#' @param value A [`character`] vector of the new cell names
#'
#' @return Updated [`PharmacoSet`]
#'
#' @importFrom CoreGx fNames<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod('fNames',
                 signature = signature(object='PharmacoSet', mDataType='character', value='character'),
                 function(object, mDataType, value) {
    callNextMethod(object=object, mDataType=mDataType, value=value)
})