#' sensitivityInfo getter
#'
#' Get the senstivity information DataFrame from a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensInf <- sensitivityInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to retrieve sensitivity experiment
#'     annotations from
#' @param dimension `character` Optional name of the dimension to extract,
#'     either 'cells' or 'drugs'. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param ... `pairlist` Additional arguments to the rowData or colData
#'     `LongTable` methods. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#'
#' @return a [`DataFrame`] with the experiment info
#'
#' @importMethodsFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the drug dose sensitivity experiment info
#'
#' @export
setMethod(sensitivityInfo, signature("PharmacoSet"),
    function(object, dimension, ...) {
        callNextMethod(object, dimension, ...)
})

#' sensitivityInfo<- setter
#'
#' Set the sensitivityInfo DataFrame in a PharmacoSet object
#'
#' @describeIn PharmacoSet Update the metadata for the treatment response
#'     experiments in the sensitivity slot.
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to update
#' @param dimension `character` Optional name of the dimension to extract,
#'     either 'cells' or 'drugs'. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param ... Additional arguments to the rowData or colData
#'     `LongTable` methods. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param value A \code{data.frame} with the new sensitivity annotations
#'
#' @return Updated \code{PharmacoSet}
#'
#' @importMethodsFrom CoreGx sensitivityInfo<-
#'
#' @export
setReplaceMethod("sensitivityInfo", signature(object="PharmacoSet", 
    value="data.frame"), function(object, dimension, ..., value) 
{
    callNextMethod(object=object, dimension=dimension, ..., value=value)    
})