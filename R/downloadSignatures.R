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
#' @param name A \code{character} string, the name of the PharmacoSet for which
#'   to download signatures. The name should match the names returned in the
#'   `PSet Name` column of `availablePSets(canonical=FALSE)`.
#' @param saveDir A \code{character} string with the folder path where the
#'   PharmacoSet should be saved. Defaults to \code{"./PSets/Sigs/"}. Will
#'   create directory if it does not exist.
library#' @param verbose \code{bool} Should status messages be printed during download.
#'   Defaults to TRUE.
#'
#' @return An array type object contaning the signatures
#'
#' @export
#' @import downloader
downloadPertSig <- function(name, saveDir=file.path(".", "PSets", "Sigs"),
    myfn=NULL, verbose=TRUE) {

    # change the download timeout since the files are big
    opts <- options()
    options(timeout=600)
    on.exit(options(opts))

    # get the annotations for available data
    pSetTable <- availablePSets(canonical=FALSE)

    # pick a signature from the list
    whichx <- match(name, pSetTable[, 3])
    if (is.na(whichx)){
        stop('Unknown Dataset. Please use the `Dataset Name` column in the
            data.frame returned by the availablePSet function to select a
            PharmacoSet')
    }
    if (!pSetTable[whichx, "type"] %in% c("perturbation", "both")){
        stop('Signatures are available only for perturbation type datasets')
    }

    if(!file.exists(saveDir)) {
        dir.create(saveDir, recursive=TRUE)
    }

    myfn <- paste(pSetTable[whichx, ]$`Dataset Name`, "_signatures.RData", sep="")

    downloader::download(paste(
        "https://www.pmgenomics.ca/bhklab/sites/default/files/downloads", myfn,
            sep='/'),
        destfile=file.path(saveDir, myfn), quiet=!verbose, mode='wb')
    sig <- load(file.path(saveDir, myfn))
    return(get(sig))
}
