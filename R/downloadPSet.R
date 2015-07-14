#' Download a PharmacoSet object
#'
#' This function allows you to download a \code{PharmacoSet} object for use with this
#' package. The \code{PharmacoSets} have been extensively curated and organised within
#' a PharacoSet class, enabling use with all the analysis tools provided in
#' \code{PharmacoGx}.
#'
#' @param name \code{Character} string, the name of the PhamracoSet to download.
#'   The available options are CGP, CCLE, and CMAP
#' @param download.method \code{Character} string, the method used by
#'   \code{download.file}. Defaults to \code{auto}.
#' @param saveDir \code{Character} string with the folder path where the
#'     PharmacoSet should be saved. Defaults to \code{'./PSets/'}. Will create
#'     directory if it does not exist.
#' @param myfn \code{character} string, the file name to save the dataset under
#' @export
#' 

downloadPSet <- function(name=c('CGP', 'CCLE', 'CMAP'), download.method = 'auto', saveDir = './PSets/', myfn=NULL) {
  
  
  name <- match.arg(name)

  if(!dir.exists(saveDir)) {
    dir.create(saveDir)
  }
  
  switch(name, 
         'CGP'={ 
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CGP.RData'
           }
         },
         'CCLE'={
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CCLE.RData'
           }
         },
         'CMAP'={
           downfn <- ''
           if (is.null(myfn)){
             myfn <- 'CMAP.RData'
           }
         }, {
           stop('Unknown Dataset. Please check the documentation of this function or the vignette for the list of available PharamcoSets.')
         })
  download.file(file.path('http://www.pmgenomics.ca/bhklab/software/pharmacogx', myfn), method=download.method, destfile = file.path(saveDir, myfn))
  pSet <- load(file.path(saveDir,myfn))
  return(pSet)
}