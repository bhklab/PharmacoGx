setOldClass('sessionInfo', sessionInfo)

#' @importFrom utils sessionInfo
.PharmacoSig <- setClass('PharmacoSig', slots=list(
            Arguments = "list",
            PSetName='character',
            DateCreated = 'character',
            SigType = 'character',
            SessionInfo = 'sessionInfo',
            Call = 'character'), contains='array')


#' Contructor for the PharmacoSig S4 class
#'
#' @param Data  of data to build the signature from
#' @param PSetName `character` vector containing name of PSet, defaults to ''
#' @param DateCreated `date` date the signature was created, defaults to `date()`
#' @param SigType `character` vector specifying whether the signature is sensitivity or perturbation, defaults to 'sensitivity'
#' @param SessionInfo `sessionInfo` object as retuned by `sesssionInfo()` function, defaults to `sessionInfo()`
#' @param Call `character` or `call` specifying the constructor call used to make the object, defaults to 'No Call Recorded'
#' @param Arguments `list` a list of additional arguments to the constructure
#'
#' @return A `PharmacoSig` object build from the provided signature data
#'
#' @examples
#' PharmacoSig()
#'
#' @export
PharmacoSig <- function(Data=array(NA, dim=c(0,0,0)), PSetName='', DateCreated=date(), SigType='sensitivity',
                        SessionInfo=sessionInfo(), Call='No Call Recorded', Arguments = list()){
    return(.PharmacoSig(Data, Arguments = Arguments, PSetName=PSetName, DateCreated=DateCreated, SigType=SigType,
                        SessionInfo=SessionInfo, Call=Call))
}


#' Show PharmacoGx Signatures
#'
#' @examples
#' data(GDSCsmall)
#' drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna",
#'              nthread=1, features = fNames(GDSCsmall, "rna")[1])
#' drug.sensitivity
#'
#' @param object \code{PharmacoSig}
#' @return Prints the PharmacoGx Signatures object to the output stream, and returns invisible NULL.
#' @export
setMethod("show", signature=signature(object='PharmacoSig'),
        function(object) {
        cat('PharmacoSet Name: ', attr(object, 'PSetName'), "\n")
        cat('Signature Type: ', attr(object, 'SigType'), "\n")
        cat("Date Created: ", attr(object, 'DateCreated'), "\n")
        cat("Number of Drugs: ", dim(object)[[2]], "\n")
        cat("Number of Genes/Probes: ", dim(object)[[1]], "\n")
           })


#' Show the Annotations of a signature object
#'
#' This funtion prints out the information about the call used to compute the drug signatures, and the session info
#' for the session in which the computation was done. Useful for determining the exact conditions used to generate signatures.
#'
#' @examples
#' data(GDSCsmall)
#' drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna",
#'              nthread=1, features = fNames(GDSCsmall, "rna")[1])
#' showSigAnnot(drug.sensitivity)
#'
#' @param object An object of the \code{PharmacoSig} Class, as
#' returned by \code{drugPerturbationSig} or \code{drugSensitivitySig}
#'
#' @return Prints the PharmacoGx Signatures annotations to the output stream, and returns invisible NULL.
#'
#' @importMethodsFrom CoreGx showSigAnnot
#' @export
setMethod("showSigAnnot", signature(object="PharmacoSig"), function(object){
  .showSigAnnotPharmacoSig(object)
})


#' @keywords internal
.showSigAnnotPharmacoSig <- function(object){
  print(object@Call)
  print(object@SessionInfo)
  return(invisible(NULL))
}
