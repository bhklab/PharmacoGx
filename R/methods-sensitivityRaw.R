#' sensitivityRaw Getter
#'
#' @describeIn PharmacoSet Retrive the raw dose and viability data from a pSet
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityRaw(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#' @return A \code{array} containing the raw sensitivity data
#'
#' @importMethodsFrom CoreGx sensitivityRaw
#' @export
setMethod("sensitivityRaw", signature("PharmacoSet"), function(object) {
    if (is(sensitivitySlot(object), 'LongTable'))
        return(.rebuildRaw(sensitivitySlot(object)))
    else
        return(callNextMethod(object))
})

#' Replicate the $raw slot in the old @sensitivity list from a LongTable
#'
#' @param longTable [`LongTable`]
#'
#' @return A 3D [`array`]
#'
#' @keywords internal
#' @noRd
.rebuildRaw <- function(longTable) {

    # Extract the information needed to reconstruct the sensitivityRaw array
    meta <- longTable$experiment_metadata[, .(rn, colKey, rowKey)]
    setkeyv(meta, c('rowKey', 'colKey'))
    dose <- assay(longTable, 'dose')
    setkeyv(dose, c('rowKey', 'colKey'))
    viab <- assay(longTable, 'viability')
    setkeyv(viab, c('rowKey', 'colKey'))

    # Join to a single wide data.table
    arrayData <- meta[dose][viab][, -c('rowKey', 'colKey')]

    # Get the data matrices
    Dose <- as.matrix(arrayData[, grep('Dose', colnames(arrayData), value=TRUE), with=FALSE])
    Viability <- as.matrix(arrayData[, grep('Viability', colnames(arrayData), value=TRUE), with=FALSE])

    array <- simplify2array(list(Dose, Viability))
    dimnames(array) <- list(arrayData$rn,
                            paste0('doses', seq_len(ncol(Dose))),
                            c('Dose', 'Viability'))
    return(array)
}

## TODO:: Make this a unit test
## TEST: all.equal(raw[dimnames(SR)[[1]],,], SR)


#' sensitivityRaw<- Setter
#'
#' @describeIn PharmacoSet Update the raw dose and viability data in a pSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityRaw(CCLEsmall) <- sensitivityRaw(CCLEsmall)
#'
#' @param object A \code{PharmacoSet} to extract the raw sensitivity data from
#' @param value A \code{array} containing the raw dose and viability data for the
#'   pSet
#'
#' @return A copy of the \code{PharmacoSet} containing the updated sensitivty data
#'
#' @importMethodsFrom CoreGx sensitivityRaw<-
#'
#' @export
setReplaceMethod('sensitivityRaw', signature("PharmacoSet", "array"), function(object, value) {
    if (is(sensitivitySlot(object), 'LongTable')) {

        ## TODO:: validate value

        longTable <- sensitivitySlot(object)

        raw <- as.data.table(value, keep.rownames=TRUE, na.rm=FALSE)

        # preprocess raw array
        ## FIXME:: refactor this into a helper, it is repeated in sensitivtySlotToLongTable
        setnames(raw, seq_len(3), c('rn', 'replicate', 'assay'))
        assayIDs <- unique(raw$assay)
        raw[, value := as.numeric(value)]
        raw[, replicate := as.integer(gsub('\\D*', '', replicate))]
        # Split value into one column for each assay (long -> wide)
        longRaw <- dcast(raw, rn + replicate ~ ..., value.var=c('value'))
        # Split assay columns into assay per replicate (wide -> wider)
        longRaw <- dcast(longRaw, rn ~ replicate, value.var=assayIDs)

        assayData <- assays(longTable, withDimnames=TRUE, key=FALSE)
        expMetadata <- assayData$experiment_metadata[, c(idCols(longTable), 'rn'), with=FALSE]
        rawData <- merge.data.table(longRaw, expMetadata, by='rn')

        doseCols <- grep('Dose_\\d+', colnames(rawData), value=TRUE)
        assayData$dose <-
            rawData[, c(idCols(longTable), doseCols), with=FALSE]
        viabilityCols <- grep('Viability_\\d+', colnames(rawData), value=TRUE)
        assayData$viability <-
            rawData[, c(idCols(longTable), viabilityCols), with=FALSE]

        assays(longTable) <- assayData
        object@sensitivity <- longTable

    } else {
        callNextMethod(object, value=value)
    }
    return(object)
})