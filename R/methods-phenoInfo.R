#' phenoInfo Getter
#'
#' Get the phenotype information for a specified molecular datatype
#'
#' @describeIn PharmacoSet Return the experiment info from the given type of
#'   molecular data in PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' phenoInf <- phenoInfo(CCLEsmall, mDataType="rna")
#'
#' @param object The [`PharmacoSet`] to retrieve rna annotations from
#' @param mDataType `character` The type of molecular data
#'
#' @return a `data.frame` with the phenotype information for the specified molecular data type
#'
#'
#' @importFrom methods callNextMethod
#' @importMethodsFrom CoreGx phenoInfo
#'
#' @export
setMethod('phenoInfo', signature='PharmacoSet', function(object, mDataType){
    callNextMethod(object=object, mDataType=mDataType)
})


#' phenoInfo<- Setter
#'
#' Update the phenotype information for a specified molecular data type in a specified pSet object
#'
#' @describeIn PharmacoSet Update the given type of molecular data experiment info in the PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' phenoInfo(CCLEsmall, mDataType='rna') <- phenoInfo(CCLEsmall, mDataType='rna')
#'
#' @param object The [`PharmacoSet`] to retrieve molecular experiment annotations from
#' @param mDataType the type of molecular data
#' @param value a `data.frame`  with the new experiment annotations
#'
#' @return The updated \code{PharmacoSet}
#'
#' @importFrom methods callNextMethod
#' @importMethodsFrom CoreGx phenoInfo<-
#'
#' @export
setReplaceMethod('phenoInfo',
                 signature = signature(object='PharmacoSet', mDataType ='character', value='DataFrame'),
                 function(object, mDataType, value){
  callNextMethod(object=object, mDataType=mDataType, value=value)
})