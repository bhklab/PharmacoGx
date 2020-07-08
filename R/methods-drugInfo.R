#' drugInfo getter
#'
#' Retrieve information from the
#'
#' @examples
#' data(CCLEsmall)
#' drugInf <- drugInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to retrieve drug info from
#'
#' @return A [`data.frame`] containg annotatations for all drugs in the object
#'
#' @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#'
#' @export
setMethod(drugInfo, signature="PharmacoSet", function(object) {
  object@drug
})

#' drugInfo Getter
#'
#' Retrieve information from the
#'
#' @examples
#' data(CCLEsmall)
#' drugInf <- drugInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to retrieve drug info from
#' @param value A [`data.frame`] containing updated drug annotations
#'
#' @return A [`PharmacoSet`] with updated drug annotations in the `@drug` slot
#'
#' @describeIn PharmacoSet Update the drug annotations
#'
#' @export
setReplaceMethod("drugInfo",
                 signature=signature(object="PharmacoSet", value="data.frame"),
                 function(object, value) {
  object@drug <- value
  object
})