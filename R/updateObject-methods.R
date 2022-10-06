#' Update the PharmacoSet class after changes in it struture or API
#'
#' @param object A `PharmacoSet` object to update the class structure for.
#'
#' @return `PharmacoSet` with update class structure.
#'
#' @md
#' @importMethodsFrom CoreGx updateObject
#' @export
setMethod("updateObject", signature("PharmacoSet"), function(object) {
    cSet <- callNextMethod(object)
    pSet <- as(cSet, "PharmacoSet")
    names(curation(pSet)) <- gsub("drug", "treatment", names(curation(pSet)))
    if ("treatment" %in% names(curation(pSet))) {
        colnames(curation(pSet)$treatment) <- gsub("treatmentid", "treatmentid",
            colnames(curation(pSet)$treatment))
    }
    validObject(pSet)
    return(pSet)
})