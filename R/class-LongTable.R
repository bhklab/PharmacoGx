#' Define an S3 Class
#'
#' Allows use of S3 methods with new S4 class. This is required to overcome
#' limitations of the `[` S4 method.
#'
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

#' Create a LongTable object from a single .csv file
#'
#' @param filePath [`character`] Path to the .csv file containing the data and
#'   metadata from which to build the `LongTable`.
#' @param colDataCols [`list`] List with two `character` vectors, the first
#'   specifying one or more columns to be used as column identifiers (e.g.,
#'   drug name columns) and the second containing any additional metadata
#'   columns related to the column identifiers.
#' @param rowDataCols [`list`] List with two `character` vectors, the first
#'   specifying one or more columns to be used as cell identifiers (e.g.,
#'   cell-line name columns) and the second containing any additional metadata
#'   columns related to the cell identifiers.
#' @param assayCols [`list`] A named list of character vectors specifying how to
#'   parse assay columns into a list of `data.table`s. Each list data.table
#'   will be named for the name of corresponding list item and contain the columns
#'   specified in the character vector of column names in each list item.
#'
#' @return A [`LongTable`] object containing one or more assays, indexed by
#'   rowID and colID.
#'
#' @import data.table
#' @export
buildLongTableFromCSV <- function(filePath, rowDataCols, colDataCols, assayCols) {

    # read in data
    tableData <- fread(filePath)

    # build drug and cell metadata tables and index by the appropriate ID
    colData <- unique(tableData[, unlist(colDataCols), with=FALSE])
    colData[, colKey := seq_len(nrow(.SD))]
    rowData <- unique(tableData[, unlist(rowDataCols), with=FALSE])
    rowData[, rowKey := seq_len(nrow(.SD))]

    # add the row and column ids to the value data
    assayData <- tableData[rowData, on=unlist(rowDataCols)][colData, on=unlist(colDataCols)]
    rm(tableData)
    assayData[, c(unlist(rowDataCols), unlist(colDataCols)) := NULL]
    setkey(assayData, rowKey, colKey)

    setkey(rowData, rowKey)
    setkey(colData, colKey)

    setnames(colData, colDataCols[[1]], paste0('drug', seq_along(colDataCols[[1]])))
    setnames(rowData, rowDataCols[[1]], paste0('cellLine', seq_along(rowDataCols[[1]])))

    # add the index columns to the different value column vectors
    # this allows the .selectDataTable helper to be more general
    .prependToVector <- function(vector, values) c(values, vector)
    assayCols <- lapply(assayCols, FUN=.prependToVector, values=c('rowKey', 'colKey'))
    assays <- lapply(assayCols, .selectDataTable, DT=assayData)

    return(LongTable(colData=colData, colIDs=paste0('drug', seq_along(colDataCols[[1]])),
                     rowData=rowData, rowIDs= paste0('cellLine', seq_along(rowDataCols[[1]])),
                     assays=assays))
}

# ---- buildLongTableFromCSV helpers

#' Select a set of column names from a data.table, returning a copy of the
#'   data.table with duplicate rows removed
#'
#' @param colNames [`character`] The column names to select from the data.table
#' @param DT [`data.table`, `data.frame`, `matrix`] An object coercible to a `data.table`.
#'   Please note rownames will be dropped by default.
#' @param keep.rownames [`logical` or `character`] Passed through to the data.table coercing if DT is not a
#'   `data.table`. If TRUE, rownames will be caputured in the `rn` column; if FALSE (default) rownames will
#'   be dropped; if `character`, rownames will be captured in a column with the same name.
#'
#' @return [`data.table`] Copy of `DT` containing only the specified columns, with duplicate rows removed.
#'
#' @import data.table
#' @keywords internal
.selectDataTable <- function(colNames, DT, keep.rownames=FALSE) {

    # validate input
    if (!is.data.table(DT)) {
        tryCatch({
            DT <- data.table(DT, keep.rownames=keep.rownames)
        }, warning=function(w) {
            warning(w)
        }, error=function(e) {
            message(e)
            stop("Argument to DT parameter must be coercible to a data.table!")
        })
    }
    if (!is.character(colnames(DT))) stop("Currently only character column ids are supported!")
    missingColumns <- setdiff(colNames, colnames(DT))
    if (length(missingColumns) > 0)
        warning(paste0("There are no columns named ", paste0(missingColumns, collapse=", "), 'in DT.
            Continuing subset without these columns.'))

    # perform subset and copy to prevent modify by refence issues
    selectedDT <- copy(unique(DT[, .SD, .SDcols=colnames(DT) %in% colNames]))

    return(selectedDT)
}

# ---- LongTable Class Methods

#' Subset method for a LongTable object.
#'
#' Allows use of the colData and rowData `data.table` objects to query based on
#'  rowID and colID, which is then used to subset all value data.tables stored
#'  in the dataList slot.
#'
#' This function is endomorphic, it always returns a LongTable object.
#'
#' @param x [`LongTable`] The object to subset.
#' @param i [`character`, `numeric`, `logical` or `expression`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all row id columns matching the vector.
#'
#'  Numeric or Logical: these select based on the rowKey from the `rowData`
#'      method for the `LongTable`.
#'
#'  Expression: Accepts valid query statements to the `data.table` i parameter,
#'      this can be used to make complex queries using the `data.table` API
#'      for the `rowData` data.table.
#'
#' @param j [`character`, `numeric`, `logical` or `expression`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all drug id columns matching the vector.
#'
#'  Numeric or Logical: these select base don the rowID from the `rowData`
#'      method for the `LongTable`.
#'
#'  Expression: Accepts valid query statements to the `data.table` i parameter,
#'      this can be used to make complex queries using the `data.table` API
#'      for the `rowData` data.table.
#'
#' @param values [`character`, `numeric` or `logical`] Optional list of value
#'      names to subset. Can be used to subset the dataList column further,
#'      returning only the selected items in the new LongTable.
#'
#' @return [`LongTable`] A new `LongTable` object subset based on the specified
#'      parameters.
#'
#' @importMethodsFrom BiocGenerics subset
#' @import data.table
#' @export
setMethod('subset', signature('LongTable'), function(x, i, j, assays) {

    longTable <- x
    rm(x)

    if (!missing(i)) {
        if (tryCatch(is.character(i), error=function(e) FALSE)) {
            select <- grep('^cellLine[:digit:]*', colnames(rowData(longTable)), value=TRUE)
            iString <- paste0(paste0(select, ' %in% ', .variableToCodeString(i)), collapse=' | ')
            i <- str2lang(iString)
        } else {
            i <- substitute(i)
        }
        rowDataSubset <- rowData(longTable)[eval(i), ]
    } else {
        rowDataSubset <- rowData(longTable)
    }

    if (!missing(j)) {
        if (tryCatch(is.character(j), error=function(e) FALSE)) {
            select <- grep('^drug[:digit:]*', colnames(colData(longTable)), value=TRUE)
            jString <- paste0(paste0(select, ' %in% ', .variableToCodeString(j)), collapse=' | ')
            j <- str2lang(jString)
        } else {
            j <- substitute(j)
        }
        colDataSubset <- colData(longTable)[eval(j), ]
    } else {
        colDataSubset <- colData(longTable)
    }

    rowKeys <- rowDataSubset$rowKey
    colKeys <- colDataSubset$colKey

    if (missing(assays)) { assays <- assayNames(longTable) }
    keepAssays <- assayNames(longTable) %in% assays

    assayData <- lapply(assays(longTable)[keepAssays],
                     FUN=.filterLongDataTable,
                     indexList=list(rowKeys, colKeys))

    return(LongTable(colData=colDataSubset, colIDs=longTable@.intern$colIDs ,
                     rowData=rowDataSubset, rowIDs=longTable@.intern$rowIDs,
                     assays=assayData, metadata=metadata(longTable)))
})


## NOTE:: Issues printing are caused by ggplot::%+% over riding crayon::%+%
#'
#'
#'
#'
#' @import crayon
#' @importFrom crayon %+%
#' @import data.table
#' @importMethodsFrom CoreGx show
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


    ## FIXME:: Why is this so slow? Also why doesn't it work?
    # --- metadata slot
    metadataLength <- length(metadata(object))
    if (metadataLength > 0) {
        metadataString <- paste0('metadata(', metadataLength, '): ')
        metadataNames <- names(metadata(object))
        metadataNamesString <-
            if (length(metadatNames) > 6)
                paste0(.collapse(head(metadataNames, 3), ' ... ', .collapse(tail(metadataNames, 3))))
            else
                .collapse(metadataNames)
        cat(yellow$bold(metadataString) %+% green(metadataNamesString), '\n')
    }
})

#' Convert an R object in a variable into a string of the code necessary to
#'   create that object
#'
#' @param variable [`Symbol`] A symbol containing an R variable
#'
#' @return [`string`] A string representation of the code necessary to
#'   reconstruct the variable.
#'
#' @keywords internal
.variableToCodeString <- function(variable) {
    codeString <- capture.output(dput(variable))
    codeString <- gsub('\"', "'", codeString)
    return(codeString)
}


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

# ---- LongTable Getter and Setter Methods

#' Retrieve the row metadata table from a LongTable object
#'
#' @param x A [`LongTable`] object to retrieve the row metadata from.
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @import data.table
#' @export
setMethod('rowData', signature(x='LongTable'), function(x) {
    return(x@rowData[, -'.rownames'])
})


#' Retrieve the column metadata table from a LongTable object
#'
#' @param x [`LongTable`]
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @import data.table
#' @export
setMethod('colData', signature(x='LongTable'), function(x) {
    return(x@colData[, -'.colnames'])
})


#' Return a list of `data.table` objects with the assay measurements from a
#'  `LongTable` object.
#'
#' @param x [`LongTable`] What to extract the assay data from.
#' @param withDimnames [`logical`] Should the returned assays be joined to
#'   the row and column identifiers (i.e., the pseudo dimnames of the object).
#'
#' @return
#'
#' @importMethodsFrom SummarizedExperiment assays
#' @export
setMethod('assays', signature(x='LongTable'), function(x, withDimnames=FALSE) {
    if (withDimnames)
        return(structure(
            lapply(assayNames(x), assay, x=x, withDimnames=withDimnames),
            .Names=assayNames(x))
            )
    else
        return(x@assays)
})

#' Retrieve an assay `data.table` object from the `assays` slot of a `LongTable`
#'    object
#'
#' @param x [`LongTable`] The `LongTable` object to get the assay from.
#' @param i [`integer`] or [`character`] vector containing the index or name
#'   of the assay, respectively.
#' @param withDimnames [`logical`] Should the dimension names be returned
#'   joined to the assay. This retrieves both the row and column identifiers
#'   and returns them attached to the
#'
#' @importMethodsFrom SummarizedExperiment assay
#' @export
setMethod('assay',
          signature(x='LongTable', i='character'),
          function(x, i, withDimnames=FALSE) {

    # validate input
    if (length(i) > 1 || !is.character(i))
        stop(paste0('Please specifying a single character assay name.'))

    keepAssay <- which(assayNames(x) == i)
    if (length(keepAssay) < 1)
        stop(paste0('There is no assay named ',
                    i,
                    ' in this LongTable. Use assayNames(longTable) for a list of
                    valid assay names.'))

    # extract the specified assay
    assayData <- assays(x)[[keepAssay]]

    # optionally join to rowData and colData
    if (withDimnames) {
        assayData <- .colIDData(x)[assayData, on='colKey'][, -'colKey']
        assayData <- .rowIDData(x)[assayData, on='rowKey'][, -'rowKey']
    }

    return(assayData)
})



# ----  LongTable Private Helpers

# ==== LongTable Accessor Methods

# ---- Private Helper Methods

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


#' Extract the row ID columns from `rowDat` of a `LongTable`
#'
#' @param object [`LongTable`] What to get the colIDData from.
#' @param key [`logical`] Should the `colKey` column also be returned? Default
#'     is TRUE.
#'
#' @export
#' @keywords internal
setMethod('.rowIDData', signature(object='LongTable'), function(object, key=TRUE) {
    colNames <- colnames(rowData(object))[.rowIDs(object)]
    keepCols <- if (key) c(colNames, 'rowKey') else colNames

    return(rowData(object)[, keepCols, with=FALSE])
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

    return(colData(object)[, keepCols, with=FALSE])
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
#' @param A [`object`] `LongTable` object to retrieve the row identifier column
#'   indexes from.
#'
#' @return [`list`] A list containing two `numeric` vectors, the first for the
#'    identifier columns for the `rowData` slot and the second with the same
#'    for the `colData` slot
#'
#' @export
#' @keywords internal
setMethod('.dimIDs', signature(object='LongTable'), function(object) {
    list(.rowIDs(object), .colIDs(object))
})

#' Retrieve the assay names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to retrieve the assay names from.
#'
#' @return [`character`] Names of the assays contained in the `LongTable`.
#'
#' @importMethodsFrom SummarizedExperiment assayNames
#' @export
setMethod('assayNames', signature(x='LongTable'), function(x) {
    return(names(assays(x)))
})

# ---- Long Table Accessor Methods

#' Get the dimensions of a `LongTable` object.
#'
#' @param x A [`LongTable`] object to retrieve dimensions for.
#'
#' @return [`numeric`] Vector of object dimensions.
#'
#' @importMethodsFrom SummarizedExperiment dim
#' @export
setMethod('dim', signature(x='LongTable'), function(x) {
    return(c(nrow(rowData(x)), nrow(colData(x))))
})

#' Get the column names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to get the column names from
#'
#' @return [`character`] Vector of column names.
#'
#' @export
setMethod('colnames', signature(x='LongTable'), function(x) {
    return(x@colData$.colnames)
})

#' Get the row names from a `LongTable` object.
#'
#' @param x A [`LongTable`] object to get the row names from
#'
#' @return [`character`] Vector of row names.
#'
#' @export
setMethod('rownames', signature(x='LongTable'), function(x) {
    return(x@rowData$.rownames)
})

#' Getter for the dimnames of a `LongTable` object
#'
#' @param x The [`LongTable`] object to retrieve the dimnames for
#'
#' @return [`list`] List with two character vectors, one for row and one for
#'     column names.
#'
#' @importMethodsFrom Biobase dimnames
#' @export
setMethod('dimnames', signature(x='LongTable'), function(x) {
    return(list(rownames(x), colnames(x)))
})

#' This method only throws an error. You may no edit the dimnames of a
#'    `LongTable` object this way.
#'
#' @param x [`LongTable`]
#'
#' @warning This function only trhows a warning then returns the original object.
#'
#' @import crayon
#' @export
setReplaceMethod('dimnames', signature(x='LongTable'), function(x, value) {
    warning(cyan$bold("The dimnames of a `LongTable` object cannot be directly
        modified. Please use the `<method_name>` instead."))
    return(x)
})

#' Getter method for the metadata slot of a `LongTable` object
#'
#' @param x The [`LongTable`] object from which to retrieve the metadata list.
#'
#' @return [`list`] The contents of the `metadata` slot of the `LongTable`
#'   object.
#'
#' @importFrom S4Vectors metadata
#' @export
setMethod('metadata', signature(x='LongTable'), function(x) {
    return(x@metadata)
})

#' Setter method for the metadata slot of a `LongTable` object
#'
#' @param x [`LongTable`] The LongTable to update
#' @param value [`list`] A list of new metadata associated with a `LongTable`
#'   object.
#'
#' @return [`LongTable`] A copy of the `LongTable` object with the `value` in
#'   the metadata slot.
#'
#' @importFrom S4Vectors `metadata<-`
#' @import crayon
#' @export
setReplaceMethod('metadata', signature(x='LongTable'), function(x, value) {
    if (!is(value, 'list'))
        stop(magenta$bold('The `metadata` slot must be a list!'))
    x@metadata <- value
    return(x)
})


#' Updates the `rowData` slot as long as the ID columns are not changed.
#'
#' @param x A [`LongTable`] object to modify.
#'
#' @return A copy of the [`LongTable`] object with the `rowData`
#'   slot updated.
#'
#' @importFrom SummarizedExperiment `rowData<-`
#' @export
#setReplaceMethod('rowData', signature(x='LongTable'), function(x, value) {
#
#    # type check input
#    if (is(value, 'data.frame'))
#        value <- data.table(value, keep.rownames=FALSE)
#    if (!is(value, 'data.table'))
#        stop(magenta$bold('Please pass a data.frame or data.table to update
#            the rowData slot. We recommend modifying the object returned by
#            rowData(x) then reassigning it with rowData(x) <- newRowData'))
#
#    # remove key column
#    if ('rowKey' %in% colnames(value)) {
#        value[, rowKey := NULL]
#        warning(cyan$bold('Dropping rowKey from replacemetn value, this
#            function will deal with mapping the rowKey automatically.'))
#    }
#
#    # assemble information to select proper update method
#    rowIDCols <- colnames(.rowIDData(x))
#    sharedRowIDCols <- intersect(rowIDCols, colnames(value))
#
#    metadataCols <- colnames(rowData(x)[, -c(rowIDCols, 'rowKey')])
#    sharedMetadataCols <- intersect(metadataCols, colnames(value))
#
#    # case where no row ids are in update
#    if (length(sharedRowIDCols) > 0) {
#        if (length(sharedMetadataCols) > 0) {
#
#        } else {
#
#        }
#    }
#
#    # case where row ids are in the updated table
#    if (all(rowIDCols %in% sharedRowIDCols)) {
#
#    } else {
#
#    }
#
#
#})



#'
#'
#' @param rowIDs [`LongTable`]
#' @param rowData [`data.table`]
#'
#' @export
#' @keywords internal
.joinOnRowIDs <- function(rowIDs, rowData) {
    if (!('rowKey' %in% colnames(rowIDs)))
        stop(magenta$bold("No rowKey column in rowIDs?"))
    row_data <- rowIDs[rowData, on=colnames(rowIDs)]
    setkeyv(row_data, 'rowKey')
    return(row_data)
}

#' Helper to determine if the rowData rowIDs are identical between a new
#'   rowData data.table object and an existing `LongTable`
#'
#' @export
#' @keywords internal
.row_IDs_are_identical <- function(rowIDs, newRowData) {
    if ('rowKey' %in% colnames(rowIDs)) rowIDs[, rowKey := NULL]
    return(all(colnames(rowIDs) %in% colnames(newRowData)) &&
        all.equal(rowIDs, newRowData[, colnames(rowIDs), with=FALSE]))
}


##'
##' @param x [`LongTable`]
##'
##' @importFrom SummarizedExperiment `colData<-`
##' @export
#setReplaceMethod('colData', signature(x='LongTable'), function(x, value) {
#
#})