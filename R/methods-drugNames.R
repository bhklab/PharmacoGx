#' drugNames Getter
#'
#' Get the names of all drugs available in a specified `PharmacoSet` object
#'
#' @describeIn PharmacoSet Return the names of the drugs used in the PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to return drug names from
#'
#' @return A `character` vector containg the names of drugs in the pSet
#'
#' @export
setMethod(drugNames, signature="PharmacoSet", function(object) {
  rownames(drugInfo(object))
})


#' drugNames<- Setter
#'
#' Set the drug names available in a PharmacoSet object
#'
#' @describeIn PharmacoSet Update the drug names used in the dataset
#'
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall) <- drugNames(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to update
#' @param value A `character` vector of the new drug names
#'
#' @return The updated [`PharmacoSet`] object
#'
#' @export
setReplaceMethod("drugNames",
                 signature=signature(object="PharmacoSet", value="character"),
                 function(object, value) {
    object <- updateDrugId(object, value)
    return(object)
})