setOldClass('sessionInfo', sessionInfo())

#' @importFrom utils sessionInfo
.PharmacoSig <- setClass('PharmacoSig', slots=list(
            
            PSetName='character',
            DateCreated = 'character',
            SigType = 'character',
            SessionInfo = 'sessionInfo',
            Call = 'character'), contains='array')

PharmacoSig <- function(Data=array(NA, dim=c(0,0,0)), PSetName='', DateCreated=date(), SigType='sensitivity', SessionInfo=sessionInfo(), Call='No Call Recorded'){

#attr(SessionInfo, 'class') <- NULL

return(.PharmacoSig(Data, PSetName=PSetName, DateCreated=DateCreated, SigType=SigType, SessionInfo=SessionInfo, Call=Call))}


#' Show PharmacoGx Signatures  
#' 
#' @param object \code{PharmacoSig}
#' 
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
#' @param Sigs An object of the \code{PharmacoSig} Class, as
#' returned by \code{drugPerturbationSig} or \code{drugSensitivitySig}
#'
#' @export
showSigAnnot <- function(Sigs){

print(Sigs@Call)
print(Sigs@SessionInfo)
} 





