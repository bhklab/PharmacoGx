###########################
## Petr Smirnov
## All rights reserved
## July 20, 2015
###########################


setOldClass('sessionInfo', sessionInfo())

.PharmacoGxSignatures <- setClass('PharmacoGxSignatures', slots=list(
			
			PSetName='character',
			DateCreated = 'character',
			SigType = 'character',
			SessionInfo = 'sessionInfo',
			Call = 'character'), contains='array')

PharmacoGxSignatures <- function(Data=array(NA, dim=c(0,0,0)), PSetName='', DateCreated=date(), SigType='sensitivity', SessionInfo=sessionInfo(), Call='No Call Recorded'){

#attr(SessionInfo, 'class') <- NULL

return(.PharmacoGxSignatures(Data, PSetName=PSetName, DateCreated=DateCreated, SigType=SigType, SessionInfo=SessionInfo, Call=Call))}


#' Show PharmacoGx Signatures  
#' 
#' @param object \code{PharmacoGxSignatures}
#' 
#' @export
setMethod("show", signature=signature(object='PharmacoGxSignatures'),
		function(object) {
		cat('PharmacoSet Name: ', attr(object, 'PSetName'), "\n")
		cat('Signature Type: ', attr(object, 'SigType'), "\n")
		cat("Date Created: ", attr(object, 'DateCreated'), "\n")
		cat("Number of Drugs: ", dim(object)[[2]], "\n")
		cat("Number of Genes/Probes: ", dim(object)[[1]], "\n")
		   })

#' Show the Annotations of a signature object
#' 
#' @param Sigs An object of the \code{PharmacoGxSignatures} Class, as
#' returned by \code{drugPerturbationSig} or \code{drugSensitivitySig}
#'
#' @export
showSigAnnot <- function(Sigs){

print(Sigs@Call)
print(Sigs@SessionInfo)
} 





