#' sensitivityProfiles PharmacoSet Method
#'
#' Get the sensitivityProfiles data.frame from a PharmacoSet object
#'
#' @describeIn PharmacoSet Return the sensitivity profile summary values for the
#'     treatment response experiment data in the sensitivity slot.
#'
#' @examples
#' data(CCLEsmall)
#' sensProf <- sensitivityProfiles(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to retrieve sensitivity experiment data from
#'
#' @return a \code{data.frame} with the experiment info
#'
#' @importFrom CoreGx sensitivityProfiles
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(sensitivityProfiles, "PharmacoSet", function(object) {
    if (is(sensitivitySlot(object), 'LongTable'))
        return(.rebuildProfiles(sensitivitySlot(object)))
    else
        return(callNextMethod(object=object))
})

#' Replicate the $profiles item in the old @sensitivity slot list object.
#'
#' @param longTable [`LongTable`]
#'
#' @export
#' @import data.table
#' @noRd
.rebuildProfiles <- function(longTable) {

    # Extract the information needed to reconstruct the sensitivityRaw array
    meta <- longTable$experiment_metadata[, .(rn, colKey, rowKey)]
    setkeyv(meta, c('rowKey', 'colKey'))

    sensProf <- assay(longTable, 'sensitivity_profiles')
    setkeyv(sensProf, c('rowKey', 'colKey'))

    profiles <- merge.data.table(meta, sensProf, all=TRUE)[, -c('colKey', 'rowKey')]
    rownames <- profiles$rn
    profiles[, rn := NULL]

    setDF(profiles, rownames=rownames)
    return(profiles)

}

## TODO:: Make this a unit test
## all.equal(prof[rownames(SP), colnames(SP)], SP)

#' sensitivityProfiles<- PharmacoSet Method
#'
#' @describeIn PharmacoSet Update the sensitivity profiles for a `PharmacoSet`
#'   object.
#'
#' @examples
#' data(GDSCsmall)
#' sensitivityProfiles(GDSCsmall) <- sensitivityProfiles(GDSCsmall)
#'
#' @param object A [`PharamcoSet`] to update.
#' @param value A [`data.frame`] with the new sensitivity profiles. If a
#'   matrix object is passed in, converted to `data.frame` before assignment.
#'
#' @return [`invisible`] Updates the `PharmacoSet` object.
#'
#' @import data.table
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature(object="PharmacoSet", value="data.frame"),
                 function(object, value) {
    if (is(sensitivitySlot(object), 'LongTable')) {
        if (!is.data.table(value)) value <- data.table(value, keep.rownames=TRUE)
        sensitivity <- sensitivitySlot(object)
        idCols <- unique(c(rowIDs(sensitivity), colIDs(sensitivity)))
        experimentMetadata <-
            assay(sensitivity, 'experiment_metadata',
                  withDimnames=TRUE, key=FALSE)[, c('rn', idCols), with=FALSE]
        setkeyv(experimentMetadata, 'rn')
        setkeyv(value, 'rn')
        value <- experimentMetadata[value][, -'rn']
        assay(object@sensitivity, 'sensitivity_profiles') <- value
    } else {
        callNextMethod(object, value=value)
    }
    return(object)
})
#' @import data.table
#' @export
setReplaceMethod("sensitivityProfiles",
                 signature(object="PharmacoSet", value="matrix"),
                 function(object, value) {
    if (is(sensitivitySlot(object), 'LongTable')) {
        if (!is.data.table(value)) value <- data.table(value, keep.rownames=TRUE)
        sensitivity <- sensitivitySlot(object)
        idCols <- unique(c(rowIDs(sensitivity), colIDs(sensitivity)))
        experimentMetadata <-
            assay(sensitivity, 'experiment_metadata',
                  withDimnames=TRUE, key=FALSE)[, c('rn', idCols), with=FALSE]
        setkeyv(experimentMetadata, 'rn')
        setkeyv(value, 'rn')
        value <- experimentMetadata[value][, -'rn']
        assay(object@sensitivity, 'sensitivity_profiles') <- value
    } else {
        callNextMethod(object, value=value)
    }
    return(object)
})

