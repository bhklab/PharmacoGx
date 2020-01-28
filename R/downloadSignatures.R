#' Download Drug Perturbation Signatures
#' 
#' This function allows you to download an array of drug perturbation
#' signatures, as would be computed by the \code{drugPerturbationSig} function,
#' for the available perturbation \code{PharmacoSets}. This function allows the
#' user to skip these very lengthy calculation steps for the datasets available,
#' and start their analysis from the already computed signatures
#' 
#' @examples
#' if (interactive()){
#' downloadPertSig("CMAP")
#' }
#'  
#' @param name \code{Character} string, the name of the PharmacoSet for which to
#'   download signatures. The name should match the names returned in the
#'   availablePSets table.
#' @param saveDir \code{Character} string with the folder path where the 
#'   PharmacoSet should be saved. Defaults to \code{"./PSets/Sigs/"}. Will
#'   create directory if it does not exist.
#' @param myfn \code{character} string, the file name to save the dataset under
#' @param verbose \code{bool} Should status messages be printed during download.
#'   Defaults to TRUE.
#' @return An array type object contaning the signatures
#' @export
#' @import downloader

downloadPertSig <- function(name, saveDir=file.path(".", "PSets", "Sigs"), myfn=NULL, verbose=TRUE) {
  
  
  pSetTable <- availablePSets()
  
  whichx <- match(name, pSetTable[,1])
  if (is.na(whichx)){
    stop('Unknown Dataset. Please use the availablePSet function for the table of available PharamcoSets.')
  }
  if (!pSetTable[whichx,"Dataset.Type"]%in%c("perturbation", "both")){
    stop('Signatures are available only for perturbation type datasets')
  }
  
  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }

  myfn <- paste(name, "_signatures.RData", sep="")
  
  downloader::download(file.path("https://www.pmgenomics.ca/bhklab/sites/default/files/downloads/", myfn), destfile=file.path(saveDir, myfn), quiet=!verbose)
  sig <- load(file.path(saveDir, myfn))
  return(get(sig))
}
