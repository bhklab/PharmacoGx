#' Download a PharmacoSet object
#'
#' This function allows you to download a \code{PharmacoSet} object for use with this
#' package. The \code{PharmacoSets} have been extensively curated and organised within
#' a PharacoSet class, enabling use with all the analysis tools provided in
#' \code{PharmacoGx}.
#'
#' @param name \code{Character} string, the name of the PharmacoSet to download.
#'   The available options are CGP, CCLE, and CMAP
#' @param saveDir \code{Character} string with the folder path where the
#'     PharmacoSet should be saved. Defaults to \code{'./PSets/'}. Will create
#'     directory if it does not exist.
#' @param myfn \code{character} string, the file name to save the dataset under
#' @param verbose \code{bool} Should status message be printed during download.
#'   Defaults to TRUE.
#' @export
#' @import downloader 

downloadPSet <- function(name=c("CGP", "CCLE", "CMAP"), saveDir=file.path(".", "PSets"), myfn=NULL, verbose=TRUE) {
  
  
  name <- match.arg(name)

  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }
  
  switch(name, 
         "CGP"={ 
             downfn <- ""
           if (is.null(myfn)){
             myfn <- "CGP.RData"
           }
         },
         "CCLE"={
             downfn <- ""
           if (is.null(myfn)){
             myfn <- "CCLE.RData"
           }
         },
         "CMAP"={
             downfn <- ""
           if (is.null(myfn)){
             myfn <- "CMAP.RData"
           }
         }, {
           stop('Unknown Dataset. Please check the documentation of this function or the vignette for the list of available PharamcoSets.')
         })
  downloader::download(file.path('http://www.pmgenomics.ca/bhklab/sites/default/files/downloads', myfn), destfile=file.path(saveDir, myfn), quiet=!verbose)
  pSet <- load(file.path(saveDir, myfn))
  return(get(pSet))
}
