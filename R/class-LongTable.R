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
                       contains='list')

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

    # handle missing parameters
    isMissing <- c(rowData=missing(rowData), rowIDs=missing(rowIDs),
                 colData=missing(colData), assays=missing(assays))

    if (any(isMissing))
        stop(.errorMsg('\nRequired parameter(s) missing: ',
            names(isMissing)[isMissing], collapse='\n\t'))

    # check parameter types and coerce or error
    if (!is(colData, 'data.table'))
        tryCatch({ setDT(colData, keep.rownames=keep.rownames) },
            error=function(e)
                stop(.errorMsg("colData must be coercible to a data.frame!")))

    if (!is(rowData, 'data.table'))
        tryCatch({ setDT(rowData, keep.rownames=keep.rownames) },
            error=function(e)
                stop(.errorMsg('rowData must be coerceible to a data.frame!')))

    isDT <- is.items(assays, FUN=is.data.table)
    isDF <- is.items(assays, FUN=is.data.frame) & !isDT
    if (!all(isDT))
        tryCatch({
            for (i in which(isDF)) setDT(assay[i], keep.rownames)
        }, error = function(e, assays) {
            message(e)
            types <- lapply(assays, typeof)
            stop(.errorMsg(
                 '\nList items are types: ',
                 types, '\nPlease ensure all items in the assays list are ',
                 'coercable to a data.frame!'), collapse=', ')
        })

    # initialize the internals object to store private metadata for a LongTable
    internals <- new.env()

    # capture row interal metadata
    if (is.numeric(rowIDs) | is.logical(rowIDs)) rowIDs <- colnames(rowData)[rowIDs]
    if (!all(rowIDs %in% colnames(rowData)))
        stop(.errorMsg('\nRow IDs not in rowData: ',
            setdiff(rowIDs, colnames(rowData)), collapse=', '))
    internals$rowIDs <- rowIDs
    lockBinding('rowIDs', internals)
    internals$rowMeta <- setdiff(colnames(rowData[, -'rowKey']), rowIDs)
    lockBinding('rowMeta', internals)

    # capture column internal metadata
    if (is.numeric(colIDs) | is.logical(colIDs))
        colIDs <- colnames(colData)[colIDs]
    if (!all(colIDs %in% colnames(colData)))
        stop(.errorMsg('\nColumn IDs not in colData: ',
            setdiff(colIDs, colnames(colData)), collapse=', '))
    internals$colIDs <- colIDs
    lockBinding('colIDs', internals)
    internals$colMeta <- setdiff(colnames(colData[, -'colKey']), colIDs)
    lockBinding('colMeta', internals)

    # Reorder columns to match the keys, this prevents issues in unit tests
    # caused by different column orders.
    setcolorder(rowData, unlist(mget(c('rowIDs', 'rowMeta'), internals)))
    setcolorder(colData, unlist(mget(c('colIDs', 'colMeta'), internals)))

    ## Assemble  the pseudo row and column names for the LongTable
    ### TODO:: Is this the slow part of the constructor?
    .pasteColons <- function(...) paste(..., collapse=':')
    rowData[, `:=`(.rownames=mapply(.pasteColons, transpose(.SD))),
        .SDcols=internals$rowIDs]
    colData[, `:=`(.colnames=mapply(.pasteColons, transpose(.SD))),
        .SDcols=internals$colIDs]

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
## FIXME:: Finish this and implement class validity methods for LongTable!
.verifyKeyIntegrity <- function(rowData, colData, assays) {
    if (!('rowKey' %in% colnames(rowData)) || !is.numeric(rowData$rowID))
        message(blue('The rowKey column is missing from rowData! Please try
            rebuilding the LongTable object with the constructor.'))
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

    # ---- class descriptions
    cat(yellow$bold$italic('< LongTable >', '\n'))
    cat(yellow$bold('dim: ', .collapse(dim(object)), '\n'))

    # --- assays slot
    assayLength <- length(assays(object))
    assaysString <- paste0('assays(', assayLength, '): ')
    assayNames <- assayNames(object)
    assayNamesString <-
        if (length(assayNames(object)) > 6)
            paste0(.collapse(head(assayNames, 3), ' ... ', .collapse(tail(assayNames, 3))))
        else
            .collapse(assayNames(object))
    cat(yellow$bold(assaysString) %+% red(assayNamesString), '\n')

    # --- rownames
    rows <- nrow(rowData(object))
    rowsString <- paste0('rownames(', rows, '): ')
    rownames <- rownames(object)
    rownamesString <-
        if (length(rownames) > 6)
            paste0(.collapse(head(rownames, 3)), ' ... ', .collapse(tail(rownames, 3)))
        else
            .collapse(rownames)
    cat(yellow$bold(rowsString) %+% green(rownamesString), '\n')

    # ---- rowData slot
    rowCols <- ncol(rowData(object))
    rowDataString <- paste0('rowData(', rowCols, '): ')
    rowColnames <- colnames(rowData(object))
    rowDataNamesString <-
        if (length(rowColnames) > 6)
            paste0(.collapse(head(rowColnames, 3)), ' ... ', .collapse(tail(rowColnames, 3)))
        else
            .collapse(rowColnames)
    cat(yellow$bold(rowDataString) %+% green(rowDataNamesString), '\n')

    # ---- colnames
    cols <- nrow(colData(object))
    colsString <- paste0('colnames(', cols, '): ')
    colnames <- colnames(object)
    colnamesString <-
        if (length(colnames) > 6)
            paste0(.collapse(head(colnames, 3)), ' ... ', .collapse(tail(colnames, 3)))
        else
            .collapse(colnames)
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


    # --- metadata slot
    metadataString <- paste0('metadata(', length(metadata(object)), '): ')
    metadataNames <- names(metadata(object))
    metadataNamesString <-
        if (length(metadataNames) > 6)
            paste0(.collapse(head(metadataNames, 3), ' ... ', .collapse(tail(metadataNames, 3))))
        else if (length(metadataNames) > 1)
            .collapse(metadataNames)
        else
            'none'
    cat(yellow$bold(metadataString) %+% green(metadataNamesString), '\n')

})


# ==== LongTable Accessor Methods

#' Get the id column names for the rowData slot of a LongTable
#'
#' @param object A [`LongTable`] to get the rowData id columns for.
#' @param data [`logical`] Should the rowData for the id columns be returned
#'     instead of the column names? Default is FALSE.
#' @param key [`logical`] Should the key column also be returned?
#'
#' @return A [`character`] vector of rowData column names if data is FALSE,
#'      otherwise a [`data.table`] with the data from the rowData id columns.
#'
#' @import data.table
#' @export
setMethod('rowIDs', signature(object='LongTable'),
    function(object, data=FALSE, key=FALSE) {

    cols <- getIntern(object, 'rowIDs')
    if (key) cols <- c(cols, 'rowKey')
    if (data) rowData(object, key=TRUE)[, ..cols] else cols
})

#' Get the id column names for the rowData slot of a LongTable
#'
#' @param object A [`LongTable`] to get the rowData metadata columns for.
#' @param data [`logical`] Should the rowData for the metadata columns be returned
#'     instead of the column names? Default is FALSE.
#' @param key [`logical`] Should the key column also be returned? Default is FALSE
#'
#' @return A [`character`] vector of rowData column names if data is FALSE,
#'      otherwise a [`data.table`] with the data from the rowData metadta columns.
#'
#' @import data.table
#' @export
setMethod('rowMeta', signature(object='LongTable'),
    function(object, data=FALSE, key=FALSE){

    cols <- getIntern(object, 'colMeta')
    if (key) cols <- c(cols, 'rowKey')
    if (data) rowData(object, key=TRUE)[, ..cols] else cols

})

#' Get the id column names for the colData slot of a LongTable
#'
#' @param object A [`LongTable`] to get the colData id columns for.
#' @param data [`logical`] Should the colData for the id columns be returned
#'     instead of the column names? Default is FALSE.
#' @param key [`logical`] Should the key column also be returned? Default is FALSE.
#'
#' @return A [`character`] vector of colData column names if data is FALSE,
#'      otherwise a [`data.table`] with the data from the colData id columns.
#'
#' @import data.table
#' @export
setMethod('colIDs', signature(object='LongTable'),
    function(object, data=FALSE, key=FALSE) {

    cols <- getIntern(object, 'colIDs')
    if (key) cols <- c(cols, 'colKey')
    if (data) colData(object, key=TRUE)[, ..cols] else cols

})

#' Get the id column names for the colData slot of a LongTable
#'
#' @param object A [`LongTable`] to get the colData metadata columns for.
#' @param data [`logical`] Should the colData for the metadata columns be returned
#'     instead of the column names? Default is FALSE.
#' @param key [`logical`] Should the key column also be returned?
#'
#' @return A [`character`] vector of colData column names if data is FALSE,
#'      otherwise a [`data.table`] with the data from the colData metadta columns.
#'
#' @import data.table
#' @export
setMethod('colMeta', signature(object='LongTable'),
    function(object, data=FALSE, key=FALSE) {

    cols <- getIntern(object, 'colMeta')
    if (key) cols <- c(cols, 'colKey')
    if (data) colData(object, key=TRUE)[, ..cols] else cols
})

#' Retrieve the value columns for the assays in a LongTable
#'
#' @param object [`LongTable`]
#' @param i Optional parameter specifying the [`character`] name or [`interger`]
#'     index of the assay to get the column names for. If missing, returns a
#'     list of value column names for all the assays.
#'
#' @return A [`list`] of `character` vectors containing the value column names for
#'     each assay if i is missing, otherwise a `character` vector of value column
#'     names for the selected assay.
#'
#' @import data.table
#' @export
setMethod('assayCols', signature(object='LongTable'),
    function(object, i) {

    colNameList <- lapply(assays(object, key=FALSE), names)
    if (!missing(i)) {
        if (length(i) > 1) stop(.errorMsg('The i parameter only accepts a ',
            'single assay name or index'))

        if ((is.numeric(i) && i < length(colNameList)) ||
            (is.character(i) && i %in% names(colNameList)))
            colNameList[[i]]
        else
            stop(.errorMsg("The specified index is invalid!"))
    } else {
        colNameList
    }

})

## TODO:: Implement a function to get the entire configuration needed to make a LongTable object