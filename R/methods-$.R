# ==== LongTable Class

#' Select an assay from a LongTable object
#'
#' @param x A [`LongTable`] object to retrieve an assay from
#' @param name [`character`] The name of the assay to get.
#'
#' @return [`data.frame`] The assay object.
#'
#' @export
setMethod('$', signature('LongTable'),
    function(x, name) {

    if (name %in% assayNames(x)) {
        x[[name]]
    } else {
        # Replicate the old list behavior of the sensitivity slot to prevent
        #   breaking any existing class methods.
        tryCatch({
            switch(name,
                   raw=.rebuildRaw(x),
                   info=.rebuildInfo(x),
                   profiles=.rebuildProfiles(x),
                   n=.rebuildN(x))
        },
        error=function(e) {
            stop(.errorMsg('There is no assay named ', name, 'in this ',
                   'LongTable, valid assays are: ', assayNames(x),
                     collapse=', '))
        })
    }
})

#' Replicate the $raw slot in the old sensitivity list
#'
#' @param longTable [`LongTable`]
#'
#' @noRd
#' @keywords internal
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

#' Replicate the $info slot in the old sensitivity list
#'
#' @param longTable [`LongTable`]
#'
#' @noRd
#' @keywords internal
.rebuildInfo <- function(longTable) {

    # Extract the information needed to reconstruct the sensitivityRaw array
    meta <- assay(longTable, 'experiment_metadata')
    setkeyv(meta, c('rowKey', 'colKey'))
    rowData <- rowData(longTable, key=TRUE)[, -'drug_cell_rep']
    setkeyv(rowData, 'rowKey')
    colData <- colData(longTable, key=TRUE)[, -'drug_cell_rep']
    setkeyv(colData, 'colKey')

    # join the tables into the original data
    info <- merge.data.table(meta, rowData, all=TRUE)
    setkeyv(info, 'colKey')
    info <- merge.data.table(info, colData, all=TRUE)[, -c('rowKey', 'colKey')]
    rownames <- info$rn
    info[, rn := NULL]

    # convert to data.frame by reference, assigning rownames
    setDF(info, rownames=rownames)

    return(info)
}

## TODO:: Make this a unit test
## all.equal(info[rownames(SI), colnames(SI)], SI

#' Replicate the $profiles slot in the old sensitivity list
#'
#' @param longTable [`LongTable`]
#'
#' @noRd
#' @keywords internal
.rebuildProfiles <- function(longTable) {

    # Extract the information needed to reconstruct the sensitivityRaw array
    meta <- longTable$experiment_metadata[, .(rn, colKey, rowKey)]
    setkeyv(meta, c('rowKey', 'colKey'))

    sensProf <- assay(longTable, 'sensitivity_profiles')
    setkeyv(sensProf, c('rowKey', 'colKey'))

    profiles <- merge.data.table(meta, sensProf, all=TRUE)
    rownames <- profiles$rn
    profiles[, rn := NULL]

    setDF(profiles, rownames=rownames)
    return(profiles)

}

## TODO:: Make this a unit test
## all.equal(prof[rownames(SP), colnames(SP)], SP)

#'
#'
#'
#'
.rebuildN <- function(longTable) {

}