#' LongTable class definition
#'
#' Define a private constructor method to be used to build a `LongTable` object.
#'
#' @param drugs [`data.table`]
#' @param cells [`data.table`]
#' @param dataList [`list`]
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
                                  .intern='environment'))

#' LongTable constructor method
#'
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
#'
#' @return [`LongTable`] object
#'
#' @import data.table
#' @export
LongTable <- function(rowData, rowIDs, colData, colIDs, assays,
                      metadata, keep.rownames=FALSE) {

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
        }, error = function(e, dataList) {
            message(e)
            types <- lapply(dataList, typeof)
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
        if (is.numeric(rowIDs) & max(rowIDs) < ncol(rowData))
            rowIDs
        else
            which(colnames(rowData) %in% rowIDs)
    lockBinding('rowIDs', internals)

    internals$colIDs <-
        if (is.numeric(colIDs) & max(colIDs) < ncol(colData))
            colIDs
        else
            which(colnames(colData) %in% colIDs)
    lockBinding('colIDs', internals)

    # Assemble the pseudo row and column names for the LongTable
    .pasteColons <- function(...) paste(..., collapse=':')
    rowData[, `:=`(.rownames=mapply(.pasteColons, transpose(.SD))), .SDcols=internals$rowIDs]
    colData[, `:=`(.colnames=mapply(.pasteColons, transpose(.SD))), .SDcols=internals$colIDs]

    return(.LongTable(rowData=rowData, colData=colData,
                      assays=assays, metadata=list(),
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
    if (!('rowID' %in% colnames(rowData)) || !is.numeric(rowData$rowID))
        stop()
    if (!('colID' %in% colnames(colData)) || !is.numeric(colData$colID))
        stop()

}

rowDataCols <- list(c("cell_line"), c("BatchID"))
colDataCols <- list(c("drugA_name", "drugB_name"))
filePath <- 'data/drug_combo_merck.csv'
assayCols <- list(dose=c('drugA Conc (µM)', 'drugB Conc (µM)'),
                  viability=c(paste0('viability', seq_len(4)),
                    'mu/muMax', 'X/X0'))

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
    colData <- unique(tableData[, .SD, .SDcols=unlist(colDataCols)])[, colKey := .I]
    rowData <- unique(tableData[, .SD, .SDcols=unlist(rowDataCols)])[, rowKey := .I]

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
            stop("Argument to DT paramater must be coercible to a data.table!")
        })
    }
    if (!is.character(colnames(DT))) stop("Currently only character column ids are supported!")
    missingColumns <- setdiff(colNames, colnames(DT))
    if (length(missingColumns) > 0)
        warning(paste0("There are no columns named ", paste0(missingColumns, collapse=", ", 'in DT.
            Continuing subset without these columns.')))

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
#' @param rowQuery [`character`, `numeric`, `logical` or `expression`]
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
#' @param columnQuery [`character`, `numeric`, `logical` or `expression`]
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
setMethod('subset', signature(x='LongTable'),
          function(x, rowQuery, columnQuery, assays) {

    longTable <- x
    rm(x)

    if (!missing(rowQuery)) {
        if (tryCatch(is.character(rowQuery), error=function(e) FALSE)) {
            ## TODO:: Add a metadata field to LongTable that stores the column and row identifier string
            select <- grep('^cellLine[:digit:]*', colnames(rowData(longTable)), value=TRUE)
            rowQueryString <- paste0(paste0(select, ' %in% ', .variableToCodeString(rowQuery)), collapse=' | ')
            rowQuery <- str2lang(rowQueryString)
        } else {
            rowQuery <- substitute(rowQuery)
        }
        rowDataSubset <- rowData(longTable)[eval(rowQuery), ]
    } else {
        rowDataSubset <- rowData(longTable)
    }

    if (!missing(columnQuery)) {
        if (tryCatch(is.character(columnQuery), error=function(e) FALSE)) {
            select <- grep('^drug[:digit:]*', colnames(colData(longTable)), value=TRUE)
            columnQueryString <- paste0(paste0(select, ' %in% ', .variableToCodeString(columnQuery)), collapse=' | ')
            columnQuery <- str2lang(columnQueryString)
        } else {
            columnQuery <- substitute(columnQuery)
        }
        colDataSubset <- colData(longTable)[eval(columnQuery), ]
    } else {
        colDataSubset <- colData(longTable)
    }

    rowIDs <- rowDataSubset$rowID
    colIDs <- colDataSubset$colID

    if (missing(assays)) { assays <- assayNames(longTable) }
    keepAssays <- assayNames(longTable) %in% assays

    assayData <- lapply(assays(longTable)[keepAssays],
                     FUN=.filterLongDataTable,
                     indexList=list(rowIDs, colIDs))

    return(LongTable(colData=colDataSubset, colIDs=longTable@.intern$colIDs ,
                     rowData=rowDataSubset, rowIDs=longTable@intern$rowIDS,
                     assays=assayData))
})

#'
#'
#'
#'
#' @importFrom crayon blue green red
#' @importMethodsFrom CoreGx show
#' @export
setMethod('show', signature(object='LongTable'), function(object) {

    ## FIXME:: Function too long. Can I refacter to a helper that prints each slot?
    .collapse <- function(...) paste0(..., collapse=' ')
    # ---- class descriptions
    cat(yellow$bold$italic('\n< LongTable >', '\n'))
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

    ## FIXME:: Why does this block add ~4s to the exectuion of this function?
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

#'
#'
#'
#'
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
    return(copy(DT[rowID %in% rowIndices & colID %in% colIndices, ]))
}

# ---- LongTable Getter and Setter Methods

#' Retrieve the row metadata table from a LongTable object
#'
#' @param x A [`LongTable`] object to retrieve the row metadata from.
#'
#' @return A [`data.table`] containing rowID, row identifiers, and row metadata.
#'
#' @importMethodsFrom SummarizedExperiment rowData
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
#' @export
setMethod('colData', signature(x='LongTable'), function(x) {
    return(x@colData[, -'.colnames'])
})


#'
#'
#'
#'
#' @importMethodsFrom SummarizedExperiment assays
#' @export
setMethod('assays', signature(x='LongTable'), function(x) {
    return(x@assays)
})

#'
#'
#'
#' @importMethodsFrom SummarizedExperiment assay
#' @export
setMethod('assay',
          signature(x='LongTable', i='character'),
          function(x, i) {

    # validate input
    if (length(i) > 1 || !is.character(i))
        stop(paste0('Please specifying a single character assay name.'))

    keepAssay <- assayNames(x) == i
    if (all(keepAssay == FALSE))
        stop(paste0('There is no assay named ',
                    i,
                    ' in this LongTable. Use assayNames(longTable) for a list of
                    valid assay names.'))

    # get specified assay
    return(assays(x)[keepAssay])
})

#'
#'
#'
#'
#' @importMethodsFrom SummarizedExperiment assayNames
#' @export
setMethod('assayNames', signature(x='LongTable'), function(x) {
    return(names(assays(x)))
})

#'
#'
#'
#' @importMethodsFrom SummarizedExperiment dim
#' @export
setMethod('dim', signature(x='LongTable'), function(x) {
    return(c(nrow(rowData(x)), nrow(colData(x))))
})

#'
#'
#'
#' @importMethodsFrom BiocGenerics colnames
setMethod('colnames', signature(x='LongTable'), function(x) {
    return(x@colData$.colnames)
})

setMethod('rownames', signature(x='LongTable'), function(x) {
    return(x@rowData$.rownames)
})

#'
#'
#' @importMethodsFrom Biobase dimnames
#' @export
setMethod('dimnames', signature(x='LongTable'), function(x) {
    return(list(rownames(x), colnames(x)))
})

#'
#'
#' @importFrom SummarizedExperiment metadata
#' @export
setMethod('metadata', signature(x='LongTable'), function(x) {
    return(x@metadata)
})