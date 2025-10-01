#' Update the PharmacoSet class after changes in it struture or API
#'
#' @param object A `PharmacoSet` object to update the class structure for.
#'
#' @return `PharmacoSet` with update class structure.
#'
#' @examples
#' data(GDSCsmall)
#' updateObject(GDSCsmall)
#'
#' @md
#' @importMethodsFrom CoreGx updateObject
#' @export
setMethod("updateObject", signature("PharmacoSet"), function(object) {
  cSet <- callNextMethod(object)
  pSet <- as(cSet, "PharmacoSet")
  names(curation(pSet)) <- gsub("drug", "treatment", names(curation(pSet)))
  if ("treatment" %in% names(curation(pSet))) {
    # Column names already match current schema; historically this was a gsub.
    # Retain the branch in case future migrations need to adjust names.
    colnames(curation(pSet)$treatment) <- colnames(curation(pSet)$treatment)
  }
  validObject(pSet)
  return(pSet)
})
