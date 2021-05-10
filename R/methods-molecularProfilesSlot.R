#' Getter for the molecular profiles slot
#'
#' @examples
#' data(CCLEsmall)
#' molProfSlot <- molecularProfilesSlot(CCLEsmall)
#'
#' @param object A [`PharmacoSet`] object to retrieve the molecular profiles slot from
#'
#' @return A `list` of `SummarizedExperiment` objects, named by molecular data type
#'
#' @describeIn PharmacoSet Getter for the molecular profiles slot
#' @importMethodsFrom CoreGx molecularProfilesSlot
#' @export
setMethod('molecularProfilesSlot', signature(object="PharmacoSet"), function(object) callNextMethod(object))

#' Setter for the molecular profiles slot
#'
#' @examples
#' data(CCLEsmall)
#' molecularProfilesSlot(CCLEsmall) <- molecularProfilesSlot(CCLEsmall)
#'
#' @param object A [`PharmacoSet`] object to update the molecular profiles slot in
#' @param value A `list` of `SummarizedExperiment` objects to update the molecular profiles slot with
#'
#' @describeIn PharmacoSet Setter for the molecular profiles slot
#'
#' @importMethodsFrom CoreGx molecularProfilesSlot<-
#' @export
setReplaceMethod('molecularProfilesSlot',
                 signature(object="PharmacoSet", value="list"),
                 function(object, value) callNextMethod(object=object, value=value)
    )