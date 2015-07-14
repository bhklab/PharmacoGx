#' Download Drug Signatures
#'
#' This function allows you to download an array of drug signatures, as would be
#' computed by the \code{drugPerturbationSig} and \code{drugSensitivitySig}
#' functions, for available perturbation and sensitivity
#' \code{PharmacoSets} respectively. This function allows the user to skip these very
#' lengthy calculation steps for the datasets available, and start their
#' analysis from the already computed signatures
#' 
#' @param name \code{Character} string, the name of the PhamracoSet for which to
#'   download signatures. The available options are CGP, CCLE, and CMAP
#' @param saveDir \code{Character} string with the folder path where the
#'     PharmacoSet should be saved. Defaults to \code{'./PSets/Sigs/'}. Will create
#'     directory if it does not exist.
#' @param myfn \code{character} string, the file name to save the dataset under
#' @export
#' @import downloader 

downloadSignatures <- function(name=c('CGP', 'CCLE', 'CMAP'), saveDir = './PSets/Sigs/', myfn=NULL) {
  
  
  name <- match.arg(name)

  if(!file.exists(saveDir)) {
    dir.create(saveDir)
  }
  
  switch(name, 
         'CGP'={ 
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CGP_signatures.RData'
           }
         },
         'CCLE'={
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CCLE_signatures.RData'
           }
         },
         'CMAP'={
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CMAP_signatures.RData'
           }
         }, {
           stop('Unknown Dataset. Please check the documentation of this function or the vignette for the list of available PharamcoSets.')
         })
  download(file.path('http://www.pmgenomics.ca/bhklab/sites/default/files/downloads/', myfn), destfile = file.path(saveDir, myfn))
  pSet <- load(file.path(saveDir,myfn))
  return(get(pSet))
}
