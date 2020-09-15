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
                       slots=list(rowData="data.table",
                                  colData="data.table",
                                  assays="list",
                                  metadata='list'))

#' LongTable constructor method
#'
#'
#'
#' @param rowData [`data.table`, `data.frame`, `matrix`] A table like object
#'   coercible to a `data.table` containing the a unique `rowID` column which
#'   is used to key assays, as well as additional row metadata to subset on.
#' @param colData [`data.table`, `data.frame`, `matrix`] A table like object
#'   coercible to a `data.table` containing the a unique `colID` column which
#'   is used to key assays, as well as additional column metadata to subset on.
#' @param assays A [`list`] containing one or more objects coercible to a
#'   `data.table`, and keyed by rowID and colID corresponding to the rowID and
#'   colID columns in colData and rowData.
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
LongTable <- function(rowData, colData, assays, keep.rownames=FALSE) {

    if (!is(colData, 'data.table')) {
        colData <- data.table(colData, keep.rownames=keep.rownames)
    }

    if (!is(rowData, 'data.table')) {
        rowData <- data.table(rowData, keep.rownames=keep.rownames)
    }

    if (!all(vapply(assays, FUN=is.data.table, FUN.VALUE=logical(1)))) {
        tryCatch({
            assays <- lapply(data, FUN=data.table, keep.rownames=keep.rownames)
        }, warning = function(w) {
            warning(w)
        }, error = function(e, dataList) {
            message(e)
            types <- lapply(dataList, typeof)
            stop(paste0('List items are types: ',
                        paste0(types, collapse=', '),
                        '\nPlease ensure all items in dataList are can be
                        coerced to data.table!'))
        })
    }

    return(.LongTable(rowData=rowData, colData=colData, assays=assays, metadata=list()))
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
#'   specifying one or more columns to be used as cell identifiers (e.g.,
#'   cell-line names columns) and the second containing any additional metadata
#'   columns related to the cell identifiers.
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
    colData <- unique(tableData[, .SD, .SDcols=unlist(colDataCols)])[, colID := .I]
    rowData <- unique(tableData[, .SD, .SDcols=unlist(rowDataCols)])[, rowID := .I]

    # add the row and column ids to the value data
    assayData <- tableData[rowData, on=unlist(rowDataCols)][colData, on=unlist(colDataCols)]
    rm(tableData)
    assayData[, c(unlist(rowDataCols), unlist(colDataCols)) := NULL]
    setkey(assayData, rowID, colID)

    setkey(rowData, rowID)
    setkey(colData, colID)

    setnames(colData, colDataCols[[1]], paste0('drug', seq_along(colDataCols[[1]])))
    setnames(rowData, rowDataCols[[1]], paste0('cellLine', seq_along(rowDataCols[[1]])))

    # add the index columns to the different value column vectors
    # this allows the .selectDataTable helper to be more general
    .prependToVector <- function(vector, values) c(values, vector)
    assayCols <- lapply(assayCols, FUN=.prependToVector, values=c('rowID', 'colID'))
    assays <- lapply(assayCols, .selectDataTable, DT=assayData)

    return(LongTable(colData=colData, rowData=rowData, assays=assays))
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

    return(LongTable(colData=colDataSubset, rowData=rowDataSubset, assays=assayData))
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
    return(x@rowData)
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
    return(x@colData)
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
#' @importMethodsFrom SummarizedExperiment dimnames
#'
setMethod('dimNames', signature(x='LongTable'), function(x) {

})
