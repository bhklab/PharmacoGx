#' molecularProfiles Getter
#'
#' Get the molecular profile data for the specified molecular data type
#'
#' @describeIn PharmacoSet Return the given type of molecular data from the PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' molProf <- molecularProfiles(CCLEsmall, "rna")
#'
#' @param object The [`PharmacoSet`] to retrieve molecular profiles from
#' @param mDataType `character` The type of molecular data
#' @param assay `character` Name of the desired assay; if excluded defaults to first assay
#'   in the SummarizedExperiment for the given mDataType
#'
#' @return a [`matrix`] of data for the given mDataType and assay
#'
#' @importMethodsFrom CoreGx molecularProfiles
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(molecularProfiles, signature="PharmacoSet", function(object, mDataType, assay){
  callNextMethod(object=object, mDataType=mDataType, assay=assay)
})


#' molecularProfiles<- Setter
#'
#' Update the molecular profile data for the specified datatype in the specified pSet object
#'
#' @describeIn PharmacoSet Update the given type of molecular data from the PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' molecularProfiles(CCLEsmall, "rna") <- molecularProfiles(CCLEsmall, "rna")
#'
#' @param object The [`PharmacoSet`] to replace molecular profiles in
#' @param mDataType The type of molecular data to be updated
#' @param assay `character` Name or index of the assay data to return
#' @param value A [`matrix`] with the new profiles
#'
#' @return Updated [`PharmacoSet`]
#'
#' @importMethodsFrom CoreGx molecularProfiles<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="PharmacoSet", mDataType ="character", assay="character", value="matrix"),
                 function(object, mDataType, assay, value){
  callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})
#' @describeIn PharmacoSet Update the given type of molecular data from the PharmacoSet
#' @export
setReplaceMethod("molecularProfiles",
                 signature = signature(object="PharmacoSet", mDataType ="character", assay="missing", value="matrix"),
                 function(object, mDataType, assay, value){
  callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})