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
#'     PharmacoSet should be saved. Defaults to \code{"./PSets/Sigs/"}. Will create
#'     directory if it does not exist.
#' @param myfn \code{character} string, the file name to save the dataset under
#' @param gene \code{bool} Should the signatures be downloaded at the gene level
#'   (TRUE) or probe level (FALSE). Defaults to TRUE.
#' @param verbose \code{bool} Should status message be printed during download.
#'   Defaults to TRUE.
#' @export
#' @import downloader 

downloadSignatures <- function(name=c("CGP", "CCLE", "CMAP"), gene=TRUE,saveDir=file.path(".", "PSets", "Sigs"), myfn=NULL, verbose=TRUE) {
  
  
  name <- match.arg(name)

  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }
  
  switch(name, 
         "CGP"={ 
           if (gene){  
             downfn <- "CGP_gene_signatures.RData"
           } else {
             downfn <- "CGP_signatures.RData"
           }
           if (is.null(myfn)){
             myfn <- downfn 
           }
         },
         "CCLE"={
           if (gene){
             downfn <- "CCLE_gene_signatures.RData"
           } else {
             downfn <- "CCLE_signatures.RData" 
           }
           if (is.null(myfn)){
             myfn <- downfn
           }
         },
         "CMAP"={
           if (gene){
             downfn <- "CMAP_gene_signatures.RData"
           } else {
             downfn <- "CMAP_signatures.RData"
           }
           if (is.null(myfn)){
             myfn <- downfn
           }
         }, {
           stop("Unknown Dataset. Please check the documentation of this function or the vignette for the list of available PharamcoSets.")
         })
  downloader::download(file.path("http://www.pmgenomics.ca/bhklab/sites/default/files/downloads/", myfn), destfile=file.path(saveDir, myfn), quiet=!verbose)
  pSet <- load(file.path(saveDir, myfn))
  return(get(pSet))
}
