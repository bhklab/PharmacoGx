#' sensitivityInfo getter
#'
#' Get the senstivity information DataFrame from a PharmacoSet object
#'
#' @examples
#' data(CCLEsmall)
#' sensInf <- sensitivityInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to retrieve sensitivity experiment
#'     annotations from
#' @param dimension [`character`] Optional name of the dimension to extract,
#'     either 'cells' or 'drugs'. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param ... [`pairlist`] Additional arguments to the rowData or colData
#'     `LongTable` methods. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#'
#' @return a [`DataFrame`] with the experiment info
#'
#' @importMethodsFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the drug dose sensitivity experiment info
#'
#' @export
setMethod(sensitivityInfo, signature("PharmacoSet"),
    function(object, dimension, ...) {
  # case where sensitivity slot is a LongTable
  if (is(sensitivitySlot(object), 'LongTable')) {
      if (!missing(dimension)) {
          switch(dimension,
              cells={ return(rowData(sensitivitySlot(object), ...)) },
              drugs={ return(colData(sensitivitySlot(object), ...)) },
              stop(.errorMsg('\n[CoreGx::sensitivityInfo] Invalid value for ',
                  'the dimension argument. Please select on of "cells" or ' ,
                  '"drugs"')))
      } else {
          return(.rebuildInfo(sensitivitySlot(object)))
      }
  # sensitivity is a list
  } else {
      if (!missing(dimension))
          warning(.warnMsg('\n[CoreGx::sensitivityInfo] The dimension argument ',
              'is only valid if the sensitivity slot contains a LongTable object.',
                  ' Ignoring the dimension and ... parameters.'))
      return(callNextMethod(object))
  }
})

#' Replicate the $info slot in the old sensitivity list from the new LongTable
#'   object
#'
#' @param longTable [`LongTable`]
#'
#' @keywords internal
#' @noRd
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



#' sensitivityInfo<- setter
#'
#' Set the sensitivityInfo DataFrame in a PharmacoSet object
#'
#' @describeIn PharmacoSet Update the metadata for the treatment response
#'     experiments in the sensitivity slot.
#'
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to update
#' @param dimension [`character`] Optional name of the dimension to extract,
#'     either 'cells' or 'drugs'. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param ... Additional arguments to the rowData or colData
#'     `LongTable` methods. Only used if the sensitivity slot contains a
#'     `LongTable` object instead of a `list`.
#' @param value A \code{data.frame} with the new sensitivity annotations
#'
#' @return Updated \code{PharmacoSet}
#'
#' @importMethodsFrom CoreGx sensitivityInfo<-
#' @import CoreGx
#' @importFrom methods callNextMethod
#'
#' @import data.table
#' @export
setReplaceMethod("sensitivityInfo",
                 signature(object="PharmacoSet", value="data.frame"),
                 function(object, dimension, ..., value) {

    if (is(sensitivitySlot(object), 'LongTable')) {
        # coerce to data.table
        if (!is.data.table(value)) value <- data.table(value, keep.rownames=TRUE)
        value[, drug_cell_rep := seq_len(.N), by=.(cellid, drugid)]
        if (missing(dimension)) {
            valueCols <- colnames(value)
            # get existing column names
            rowDataCols <- colnames(rowData(object@sensitivity))
            colDataCols <- colnames(colData(object@sensitivity))
            experimentCols <- colnames(assay(object@sensitivity,
                'experiment_metadata', withDimnames=TRUE, key=FALSE))
            # drop any value columns that don't already exist
            hasValueCols <- valueCols %in% c(rowDataCols, colDataCols, experimentCols)
            if (!all(hasValueCols))
                warning(.warnMsg('\n[PharmacoGx::sensitivityInfo<-] Dropping ',
                    'columns ', .collapse(valueCols[!hasValueCols]), ' from ',
                    'value. Currently this setter only allows modifying ',
                    'existing columns when @sensitivity is a LongTable. For ',
                    'more fine grained updates please use the dimension ',
                    'argument.',  collapse=', '))
            # update the object
            message("Doing rowData")
            rowData(object@sensitivity, ...) <-
                unique(value[, .SD, .SDcols=valueCols %in% rowDataCols])
            message("Doing colData")
            colData(object@sensitivity, ...) <-
                unique(value[, .SD, .SDcols=valueCols %in% colDataCols])
            message('Doing experment_metadata')
            assay(object@sensitivity, 'experiment_metadata') <-
                unique(value[, .SD, .SDcols=valueCols %in% experimentCols])
        } else {
            switch(dimension,
                cells={ rowData(object@sensitivity, ...) <- value },
                drugs={ colData(object@sensitivity, ...) <- value },
                experiments={ assay(object@sensitivity, 'experiment_metadata') <- value },
                stop(.errorMsg('\n[PharmacoGx::sensitivityInfo<-] Invalid argument to',
                    'dimension parameter. Please choose one of "cells" or "drugsA"')))
        }
    } else {
        if (!missing(dimension))
            warning(.warnMsg('\n[PharmacoGx::sensitivityInfo<-] The dimension argument ',
                'is only valid if the sensitivity slot contains a LongTable object. ',
                'Ignoring dimension and ... parameters.'))
        object <- callNextMethod(object, value=value)
    }
    return(object)
})