#' S3 Class for LongTable S4 Class
#'
#' Allows use of S3 methods with new S4 class. This is required to overcome
#' limitations of the `[` S4 method.
#'
#' @export
setOldClass('long.table')

#' LongTable class definition
#'
#' Define a private constructor method to be used to build a `LongTable` object.
#'
#' @param drugs [`data.table`]
#' @param cells [`data.table`]
#' @param assays [`list`]
#' @param metadata [`list`]
#'
#' @return [`LongTable`] object containing the assay data from a
#'
#' @import data.table
#' @keywords internal
.LongTable <- setClass("LongTable",
                       slots=list(rowData='data.table',
                                  colData='data.table',
                                  assays='list',
                                  metadata='list',
                                  .intern='environment'),
                       contains='long.table')

#' LongTable constructor method
#'
#' Constrcuts a long table #FIXME:: Better description
#'
#' @param rowData [`data.table`, `data.frame`, `matrix`] A table like object
#'   coercible to a `data.table` containing the a unique `rowID` column which
#'   is used to key assays, as well as additional row metadata to subset on.
#' @param rowIDs [`character`, `integer`] A vector specifying
#'   the names or integer indexes of the row data identifier columns. These
#'   columns will be pasted together to make up the row.names of the
#'   `LongTable` object.
#' @param colData [`data.table`, `data.frame`, `matrix`] A table like object
#'   coercible to a `data.table` containing the a unique `colID` column which
#'   is used to key assays, as well as additional column metadata to subset on.
#' @param colIDs [`character`, `integer`] A vector specifying
#'   the names or integer indexes of the col data identifier columns. These
#'   columns will be pasted together to make up the col.names of the
#'   `LongTable` object.
#' @param assays A [`list`] containing one or more objects coercible to a
#'   `data.table`, and keyed by rowID and colID corresponding to the rowID and
#'   colID columns in colData and rowData.
#' @param metadata A [`list`] of metadata associated with the `LongTable`
#'   object being constructed
#' @param keep.rownames [`logical` or `character`] Logical: whether rownames
#'   should be added as a column if coercing to a `data.table`, default is FALSE.
#'   If TRUE, rownames are added to the column 'rn'. Character: specify a custom
#'   column name to store the rownames in.
#'
#' @return [`LongTable`] object
#'
#' @import data.table
#' @export
LongTable <- function(rowData, rowIDs, colData, colIDs, assays,
                      metadata=list(), keep.rownames=FALSE) {

    ## TODO:: Handle missing parameters

    if (!is(colData, 'data.table')) {
        colData <- data.table(colData, keep.rownames=keep.rownames)
    }

    if (!is(rowData, 'data.table')) {
        rowData <- data.table(rowData, keep.rownames=keep.rownames)
    }

    if (!all(vapply(assays, FUN=is.data.table, FUN.VALUE=logical(1)))) {
        tryCatch({
            assays <- lapply(assays, FUN=data.table, keep.rownames=keep.rownames)
        }, warning = function(w) {
            warning(w)
        }, error = function(e, assays) {
            message(e)
            types <- lapply(assays, typeof)
            stop(paste0('List items are types: ',
                        paste0(types, collapse=', '),
                        '\nPlease ensure all items in the assays list are
                        coerced to data.tables!'))
        })
    }

    # Initialize the .internals object to store private metadata for a LongTable
    internals <- new.env()

    ## TODO:: Implement error handling
    internals$rowIDs <-
        if (is.numeric(rowIDs) && max(rowIDs) < ncol(rowData))
            rowIDs
        else
            which(colnames(rowData) %in% rowIDs)
    lockBinding('rowIDs', internals)

    internals$colIDs <-
        if (is.numeric(colIDs) && max(colIDs) < ncol(colData))
            colIDs
        else
            which(colnames(colData) %in% colIDs)
    lockBinding('colIDs', internals)

    # Assemble the pseudo row and column names for the LongTable
    .pasteColons <- function(...) paste(..., collapse=':')
    rowData[, `:=`(.rownames=mapply(.pasteColons, transpose(.SD))), .SDcols=internals$rowIDs]
    colData[, `:=`(.colnames=mapply(.pasteColons, transpose(.SD))), .SDcols=internals$colIDs]

    return(.LongTable(rowData=rowData, colData=colData,
                      assays=assays, metadata=metadata,
                      .intern=internals))
}

#' Ensure that all rowID and colID keys are valid
#'
#' @param rowData [`data.table`]
#' @param colData [`data.table`]
#' @param assays [`list`]
#'
#' @keywords internal
.verifyKeyIntegrity <- function(rowData, colData, assays) {
    if (!('rowKey' %in% colnames(rowData)) || !is.numeric(rowData$rowID))
        message(blue('The rowKey column is missing from rowData! Please try rebuilding the LongTable object with the constructor.'))
    if (!('colKey' %in% colnames(colData)) || !is.numeric(colData$colID))
        stop()

}

# ---- LongTable Class Methods

## NOTE:: Issues printing are caused by ggplot::%+% over riding crayon::%+%
#' Show method for the LongTable class
#'
#' @param object A [`LongTable`] object to print the results for.
#'
#' @importFrom crayon %+% yellow red green blue cyan magenta
#' @import data.table
#' @export
setMethod('show', signature(object='LongTable'), function(object) {

    ## FIXME:: Function too long. Can I refacter to a helper that prints each slot?
    .collapse <- function(...) paste0(..., collapse=' ')

    #`%+%` <- crayon::`%+%`

    # ---- class descriptions
    cat(yellow$bold$italic('< LongTable >', '\n'))
    cat(yellow$bold('dim: ', .collapse(dim(object)), '\n'))

    # --- assays slot
    assayLength <- length(assays(object))
    assaysString <- paste0('assays(', assayLength, '): ')
    assayNames <- assayNames(object)
    assayNamesString <-
        if (length(assayNames(object)) > 6) {
            paste0(.collapse(head(assayNames, 3), ' ... ', .collapse(tail(assayNames, 3))))
        } else {
            .collapse(assayNames(object))
        }
    cat(yellow$bold(assaysString) %+% red(assayNamesString), '\n')

    # --- rownames
    rows <- nrow(rowData(object))
    rowsString <- paste0('rownames(', rows, '): ')
    rownames <- rownames(object)
    rownamesString <-
        if (length(rownames) > 6) {
            paste0(.collapse(head(rownames, 3)), ' ... ', .collapse(tail(rownames, 3)))
        } else {
            .collapse(rownames)
        }
    cat(yellow$bold(rowsString) %+% green(rownamesString), '\n')

    # ---- rowData slot
    rowCols <- ncol(rowData(object))
    rowDataString <- paste0('rowData(', rowCols, '): ')
    rowColnames <- colnames(rowData(object))
    rowDataNamesString <-
        if (length(rowColnames) > 6) {
            paste0(.collapse(head(rowColnames, 3)), ' ... ', .collapse(tail(rowColnames, 3)))
        } else {
            .collapse(rowColnames)
        }
    cat(yellow$bold(rowDataString) %+% green(rowDataNamesString), '\n')

    # ---- colnames
    cols <- nrow(colData(object))
    colsString <- paste0('colnames(', cols, '): ')
    colnames <- colnames(object)
    colnamesString <-
        if (length(colnames) > 6) {
            paste0(.collapse(head(colnames, 3)), ' ... ', .collapse(tail(colnames, 3)))
        } else {
            .collapse(colnames)
        }
    cat(yellow$bold(colsString) %+% green(colnamesString), '\n')

    # ---- colData slot
    colCols <- ncol(colData(object))
    colDataString <- paste0('colData(', colCols, '): ')
    colColnames <- colnames(colData(object))
    colDataNamesString <-
        if (length(colColnames) > 6)
            paste0(.collapse(head(colColnames, 3)), ' ... ', .collapse(tail(colColnames, 3)))
        else
            .collapse(colColnames)
    cat(yellow$bold(colDataString) %+% green(colDataNamesString), '\n')


    ## FIXME:: Why is this so slow? Also why doesn't it work?
    # --- metadata slot
    metadataString <- paste0('metadata(', length(metadata(object)), '): ')
    metadataNames <- names(metadata(object))
    metadataNamesString <-
        if (length(metadataNames) > 6) {
            paste0(.collapse(head(metadataNames, 3), ' ... ', .collapse(tail(metadataNames, 3))))
        }
        else if (length(metadataNames) > 1) {
            .collapse(metadataNames)
        }
        else {
            'none'
        }
    cat(yellow$bold(metadataString) %+% green(metadataNamesString), '\n')

})

#' Filter a data.table object based on the rowID and colID columns
#'
#' @param DT [`data.table`] Object with the columns rowID and colID, preferably
#'  as the key columns.
#' @param indexList [`list`] Two integer vectors, one indicating the rowIDs and
#'  one indicating the colIDs to filter the `data.table` on.
#'
#' @return [`data.table`] A copy of `DT` subset on the row and column IDs specified
#'  in `indexList`.
#'
#' @import data.table
#' @keywords internal
.filterLongDataTable <- function(DT, indexList) {

    # validate input
    if (length(indexList) > 2)
        stop("This object is 2D, please only pass in two ID vectors, one for
             rows and one for columns!")

    if (!all(vapply(unlist(indexList), is.numeric, FUN.VALUE=logical(1))))
        stop('Please ensure indexList only contains integer vectors!')

    # extract indices
    rowIndices <- indexList[[1]]
    colIndices <- indexList[[2]]

    # return filtered data.table
    return(copy(DT[rowKey %in% rowIndices & colKey %in% colIndices, ]))
}

# ==== LongTable Accessor Methods

# ---- Private Helper Methods

# generics

#' Return the identifiers for the column meta data in an object
#'
#' @export
setGeneric('.colIDData', function(object, ...) standardGeneric('.colIDData'))

#' Return the identifiers for the row metadata columns in an object
#'
#'
#' @export
setGeneric('.rowIDData', function(object, ...) standardGeneric('.rowIDData'))


#' Private method to retrieve the .colIDs property from an object
#'
#' @export
setGeneric('.colIDs', function(object, ...) standardGeneric('.colIDs'))

#' Private method to retrieve the .rowIDs property from an object
#'
#' @export
#' @keywords internal
setGeneric('.rowIDs', function(object, ...) standardGeneric('.rowIDs'))

#' Private method to retrieve the both the .rowIDs and .colIDs properties from
#'   an object in a list.
#'
#' @param object [`any`] An object with a .intern slot with the items .rowIDs
#'   and .colIDs
#'
#' @return A [`list`] v
#'
#' @export
#' @keywords internal
setGeneric('.dimIDs', function(object, ...) standardGeneric('.dimIDs'))

# methods

#' Extract the row ID columns from `rowDat` of a `LongTable`
#'
#' @param object [`LongTable`] What to get the colIDData from.
#' @param key [`logical`] Should the `colKey` column also be returned? Default
#'     is TRUE.
#'
#' @export
#' @keywords internal
setMethod('.rowIDData', signature(object='LongTable'), function(object, key=TRUE) {
    colNames <- colnames(rowData(object, key))[.rowIDs(object)]
    keepCols <- if (key) c(colNames, 'rowKey') else colNames

    return(rowData(object, key)[, keepCols, with=FALSE])
})

#' Extract the column ID columns from colData of a LongTable
#'
#' @param object [`LongTable`] What to get the colIDData from
#' @param key [`logical`] Should the colKey also be returned? Default is TRUE.
#'
#' @export
#' @keywords internal
setMethod('.colIDData', signature(object='LongTable'), function(object, key=TRUE) {
    colNames <- colnames(colData(object))[.colIDs(object)]
    keepCols <- if (key) c(colNames, 'colKey') else colNames

    return(colData(object, key)[, keepCols, with=FALSE])
})

#' Developer accessor method to determine which columns hold the column identifiers
#'   in the `colData` slot of a `LongTable` object.
#'
#' @param A [`object`] `LongTable` object to retrieve the column identifier
#'    column indexes from.
#'
#' @return [`numeric`] A numeric vector of integer column indexes for the ID
#'   columns in `colData`.
#'
#' @export
#' @keywords internal
setMethod('.colIDs', signature(object='LongTable'), function(object) {
    object@.intern$colIDs
})

#' Developer accessor method to determine which columns hold the row identifiers
#'   in the `rowData` slot of a `LongTable` object.
#'
#' @param A [`object`] `LongTable` object to retrieve the row identifier column
#'   indexes from.
#'
#' @return [`numeric`] A numeric vector of integer column indexes for the ID
#'   columns in `rowData`.
#'
#' @export
#' @keywords internal
setMethod('.rowIDs', signature(object='LongTable'), function(object) {
    object@.intern$rowIDs
})

#' Developer accessor method to determine which columns hold the row identifiers
#'   in the `rowData` slot of a `LongTable` object.
#'
#' @param A [`LongTable`] object to retrieve the row identifier column
#'   indexes from.
#'
#' @return A [`list`] containing two `numeric` vectors, the first for the
#'    identifier columns for the `rowData` slot and the second with the same
#'    for the `colData` slot
#'
#' @export
#' @keywords internal
setMethod('.dimIDs', signature(object='LongTable'), function(object) {
    list(.rowIDs(object), .colIDs(object))
})