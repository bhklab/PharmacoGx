#' Download Drug Perturbation Signatures
#' 
#' This function allows you to download an array of drug perturbation
#' signatures, as would be computed by the `drugPerturbationSig` function,
#' for the available perturbation `PharmacoSets`. This function allows the
#' user to skip these very lengthy calculation steps for the datasets available,
#' and start their analysis from the already computed signatures
#' 
#' @examples
#'
#' \dontrun{
#'     if (interactive()) downloadPertSig("CMAP_2016")
#' }
#' 
#' @param name A `character(1)` string, the name of the PharmacoSet for which
#'   to download signatures. The name should match the names returned in the
#'   `PSet Name` column of `availablePSets(canonical=FALSE)`.
#' @param saveDir A `character(1)` string with the folder path where the
#'   PharmacoSet should be saved. Defaults to `"./PSets/Sigs/"`. Will
#'   create directory if it does not exist.
#' @param fileName `character(1)` What to name the downloaded file. Defaults
#' to '`name`_signature.RData' when excluded.
#' @param verbose `logical(1)` Should `downloader` show detailed messages?
#' @param ... `pairlist` Force subsequent arguments to be named.
#' @param myfn `character(1)` A deprecated version of `fileName`. Still works
#' for now, but will be deprecated in future releases.
#'
#' @return An array type object contaning the signatures
#'
#' @export
#' @importFrom CoreGx .warning .funContext
#' @import downloader
downloadPertSig <- function(name, saveDir=file.path(".", "PSets", "Sigs"), 
    fileName, verbose=TRUE, ..., myfn) 
{
    funContext <- .funContext('::downloadPertSig')
    if (missing(fileName) && !missing(myfn)) {
        .warning(funContext, 'The `myfn` parameter is being deprecated in 
            favour of `fileName`. It still works for now, but will be retired
            in a future release.')
        fileName <- myfn
    }

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

    if (missing(fileName)) {
        fileName <- paste(pSetTable[whichx, ]$`Dataset Name`, 
            "_signatures.RData", sep="")
    }

    downloader::download(
        paste("https://www.pmgenomics.ca/bhklab/sites/default/files/downloads", 
            fileName, sep='/'), 
        destfile=file.path(saveDir, fileName), quiet=!verbose, mode='wb')

    sig <- load(file.path(saveDir, fileName))

    return(get(sig))
}
