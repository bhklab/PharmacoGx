#' Return a table of PharmacoSets available for download
#' 
#' The function fetches a table of all PharmacoSets available for download from the PharmacoGx server. 
#' The table includes the names of the PharamcoSet, the types of data available
#' in the object, and the date of last update.
#'
#' @examples
#' if (interactive()){
#' availablePSets()
#' }
#' 
#' @param saveDir \code{character} Directory to save the table of PSets
#' @param myfn \code{character} The filename for the table of PSets
#' @param verbose \code{bool} Should status messages be printed during download.
#' @return A data.frame with details about the available PharmacoSet objects
#' @export
#' @import downloader
#' @importFrom utils read.table write.table
availablePSets <- function(saveDir=file.path(".", "PSets"), myfn="PSets.csv", verbose=TRUE){
  
  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }  
  
  downloader::download("https://www.pmgenomics.ca/bhklab/sites/default/files/downloads/PSets.csv", destfile=file.path(saveDir, myfn), quiet=!verbose)
  pSetTable <- read.table(file.path(saveDir, myfn))
  return(pSetTable)
}


#' Download a PharmacoSet object
#'
#' This function allows you to download a \code{PharmacoSet} object for use with this
#' package. The \code{PharmacoSets} have been extensively curated and organised within
#' a PharacoSet class, enabling use with all the analysis tools provided in
#' \code{PharmacoGx}.
#' 
#' @examples
#' if (interactive()){
#' downloadPSet("CMAP")
#' }
#' 
#' @param name \code{Character} string, the name of the PhamracoSet to download.
#' @param saveDir \code{Character} string with the folder path where the
#'     PharmacoSet should be saved. Defaults to \code{'./PSets/'}. Will create
#'     directory if it does not exist.
#' @param pSetFileName \code{character} string, the file name to save the dataset under
#' @param verbose \code{bool} Should status messages be printed during download.
#'   Defaults to TRUE.
#' @return A PSet object with the dataset, downloaded from our server
#' @export
#' @import downloader 

downloadPSet <- function(name, saveDir=file.path(".", "PSets"), pSetFileName=NULL, verbose=TRUE) {
  
  pSetTable <- availablePSets(saveDir=saveDir)
  
  whichx <- match(name, pSetTable[,1])
  if (is.na(whichx)){
    stop('Unknown Dataset. Please use the availablePSets() function for the table of available PharamcoSets.')
  }
  
  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }  
  
  if(is.null(pSetFileName)){
    pSetFileName <- paste(pSetTable[whichx,"PSet.Name"], ".RData", sep="")
  }
  if(!file.exists(file.path(saveDir, pSetFileName))){
    downloader::download(url = as.character(pSetTable[whichx,"URL"]), destfile=file.path(saveDir, pSetFileName), quiet=!verbose)
  }
  pSet <- load(file.path(saveDir, pSetFileName))
  return(get(pSet))
}

#' @importFrom utils read.table write.table
.createPSetEntry <- function(pSet, outfn) {
  
  if(file.exists(outfn)){
    pSetTable <- read.table(outfn, as.is=TRUE)
    newrow <- c(pSetName(pSet), pSet@datasetType, paste(names(pSet@molecularProfiles), collapse="/"), pSet@annotation$dateCreated, NA)
    pSetTable <- rbind(pSetTable, newrow)
    rownames(pSetTable) <- pSetTable[,1]
    write.table(pSetTable, file=outfn)
  } else {
    newrow <- c(pSetName(pSet), pSet@datasetType, paste(names(pSet@molecularProfiles), collapse="/"), pSet@annotation$dateCreated, NA)
    pSetTable <- t(matrix(newrow))
    colnames(pSetTable) <- c("PSet.Name","Dataset.Type","Available.Molecular.Profiles","Date.Updated","URL")
    rownames(pSetTable) <- pSetTable[,1]
    write.table(pSetTable, file=outfn)
  }
}
