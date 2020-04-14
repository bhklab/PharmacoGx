#' Return a table of PharmacoSets available for download
#' 
#' The function fetches a table of all PharmacoSets available for download. 
#' The table includes the dataset names, version information for the data in the PSet,
#' the date of last update, the name of the PSet, and references for the data contained within,
#' a DOI for the data, and a direct download link. Download can also be done using the downloadPSet
#' function.  
#'
#' @examples
#' if (interactive()){
#' availablePSets()
#' }
#' 
#' @return A data.frame with details about the available PharmacoSet objects
#' @export
#' @import jsonlite
availablePSets <- function(){
  
  avail.psets <- fromJSON("http://orcestra.ca/api/psets/available")

  pSetTable <- data.frame("Dataset Name" = avail.psets$dataset$name,
                          "Date Created" = avail.psets$dateCreated,
                          "PSet Name" = avail.psets$name,
                          avail.psets$dataset$versionInfo, 
                          "DOI" = avail.psets$doi,
                          "Download" = avail.psets$downloadLink, stringsAsFactors = FALSE, check.names = FALSE)

  return(pSetTable)
}


#' Download a PharmacoSet object
#'
#' This function allows you to download a \code{PharmacoSet} object for use with this
#' package. The \code{PharmacoSets} have been extensively curated and organised within
#' a PharacoSet class, enabling use with all the analysis tools provided in
#' \code{PharmacoGx}. User \code{availablePSets} to discover which PSets are available.
#' 
#' @examples
#' if (interactive()){
#' downloadPSet("CTRPv2")
#' }
#' 
#' @param name \code{Character} string, the name of the PhamracoSet to download. 
#' Note that this is not the dataset name, but the PSet name - dataset names are not guaranteed
#' to be unique. 
#' @param saveDir \code{Character} string with the folder path where the
#'     PharmacoSet should be saved. Defaults to \code{'./PSets/'}. Will create
#'     directory if it does not exist.
#' @param pSetFileName \code{character} string, the file name to save the dataset under
#' @param verbose \code{bool} Should status messages be printed during download.
#'   Defaults to TRUE.
#' @return A PSet object with the dataset
#' @export
#' @import downloader 

downloadPSet <- function(name, saveDir=file.path(".", "PSets"), pSetFileName=NULL, verbose=TRUE) {
  
  pSetTable <- availablePSets()
  
  whichx <- match(name, pSetTable[,"PSet Name"])
  if (is.na(whichx)){
    stop('Unknown Dataset. Please use the availablePSets() function for the table of available PharamcoSets.')
  }
  
  if(!file.exists(saveDir)) {
    dir.create(saveDir, recursive=TRUE)
  }  
  
  if(is.null(pSetFileName)){
    pSetFileName <- paste(pSetTable[whichx,"PSet Name"], ".rds", sep="")
  }
  if(!file.exists(file.path(saveDir, pSetFileName))){
    downloader::download(url = as.character(pSetTable[whichx,"Download"]), destfile=file.path(saveDir, pSetFileName), quiet=!verbose)
  }
  pSet <- readRDS(file.path(saveDir, pSetFileName))
  return(pSet)
}

#' @importFrom utils read.table write.table
.createPSetEntry <- function(pSet, outfn) {
  
  if(file.exists(outfn)){
    pSetTable <- read.table(outfn, as.is=TRUE)
    newrow <- c(name(pSet), pSet@datasetType, paste(names(pSet@molecularProfiles), collapse="/"), pSet@annotation$dateCreated, NA)
    pSetTable <- rbind(pSetTable, newrow)
    rownames(pSetTable) <- pSetTable[,1]
    write.table(pSetTable, file=outfn)
  } else {
    newrow <- c(name(pSet), pSet@datasetType, paste(names(pSet@molecularProfiles), collapse="/"), pSet@annotation$dateCreated, NA)
    pSetTable <- t(matrix(newrow))
    colnames(pSetTable) <- c("PSet.Name","Dataset.Type","Available.Molecular.Profiles","Date.Updated","URL")
    rownames(pSetTable) <- pSetTable[,1]
    write.table(pSetTable, file=outfn)
  }
}
