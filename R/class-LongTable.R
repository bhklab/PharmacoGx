#' LongTable class for storing
#'
#'
#' @param drugs [`data.table`]
#' @param cells [`data.table`]
#' @param dataList [`list`]
#'
#' @return [`LongTable`] object containing the viability data from a treatment response experiment, indexes by a primary
#'  key
#'
#' @keywords internal
.LongTable <- setClass("LongTable",
                       slots=list(drugs="data.table",
                                  cells="data.table",
                                  dataList="list"))

#'
#'
#'
#'
#'
#'
LongTable <- function(drugs, cells, dataList) {

    if (!is(drugs, 'data.table')) {
        drugs <- data.table(drugs, keep.rownames='rownames')
    }

    if (!is(cells, 'data.table')) {
        cells <- data.table(cells, keep.rownames='rownames')
    }

    if (!all(vapply(dataList, FUN=is.data.table, FUN.VALUE=logical(1)))) {
        tryCatch({
            dataList <- lapply(data, FUN=data.table, keep.rownames='rownames')
        }, warning = function(w) {
            warning(w)
        }, error = function(e, dataList) {
            message(e)
            types <- lapply(dataList, typeof)
            stop(paste0("List items are types: ",
                        paste0(types, collapse=', '),
                        '\nPlease ensure all items in dataList are can be coerced to data.table!'))
        })
    }

    return(.LongTable(drugs=drugs, cells=cells, dataList=dataList))
}

cellCols <- list(c("cell_line"), c("BatchID"))
drugCols <- list(c("drugA_name", "drugB_name"))
filePath <- 'data/drug_combo_merck.csv'
valueCols <- list(dose=c('drugA Conc (µM)', 'drugB Conc (µM)'),
                  viability=c(paste0('viability', seq_len(4)), 'mu/muMax', 'X/X0'))

#' Create a LongTable object from a single .csv file
#'
#' @param filePath [`character`] Path to the .csv file containing the data and metadata from which to build the
#'  LongTable.
#' @param cellCols [`list`] List with two `character` vectors, the first specifying one or more columns to be used as
#'   cell identifiers (e.g., cell-line names columns) and the second containing any additional metadata columns related
#'   to the cell identifiers.
#' @param drugCols [`list`] List with two `character` vectors, the first specifying one or more columns to be used as
#'   cell identifiers (e.g., cell-line names columns) and the second containing any additional metadata columns related
#'   to the cell identifiers.
#' @param valueCols [`list`] A named list of character vectors specifying how to parse value columns into a list of
#'   `data.table`s. Each list data.table will be named for the name of corresponding list item and contain the columns
#'   specified in the character vector of column names in each list item.
#'
#' @return [`LongTable`] A long table object containing
#'
#' @import data.table
#' @export
buildLongTableFromCSV <- function(filePath, cellCols, drugCols, valueCols) {

    # read in data
    rawData <- fread(filePath)

    # build drug and cell metadata tables and index by the appropriate ID
    drugData <- unique(rawData[, .SD, .SDcols=unlist(drugCols)])[, colID := .I]
    cellData <- unique(rawData[, .SD, .SDcols=unlist(cellCols)])[, rowID := .I]

    # add the row and column ids to the value data
    valueData <- rawData[drugData, on=unlist(drugCols)][cellData, on=unlist(cellCols)]
    rm(rawData)
    valueData[, c(unlist(cellCols), unlist(drugCols)) := NULL]
    setkey(valueData, rowID, colID)
    setkey(drugData, colID)
    setkey(cellData, rowID)

    setnames(drugData, drugCols[[1]], paste0('drug', seq_along(drugCols[[1]])))
    setnames(cellData, cellCols[[1]], paste0('cellLine', seq_along(cellCols[[1]])))

    # add the index columns to the different value column vectors
    # this allows the .selectDataTable helper to be more general
    .prependToVector <- function(vector, values) c(values, vector)
    valueCols <- lapply(valueCols, FUN=.prependToVector, values=c('rowID', 'colID'))
    dataList <- lapply(valueCols, .selectDataTable, DT=valueData)

    return(LongTable(drugs=drugData, cells=cellData, dataList=dataList))
}

# ---- buildLongTableFromCSV helpers

#' Select a set of column names from a data.table, returning a copy of the data.table with duplicate rows removed
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
        }, warning=function(w){
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
          function(x, rowQuery, columnQuery, values) {

    longTable <- x
    rm(x)

    if (!missing(rowQuery)) {
        #if (tryCatch(is.character(rowQuery), error=function(e) FALSE)) {
        #
        #    select <- grep('^drug[\s')
        #    cellDTSubset <-
        #} else {
            rowQuery <- substitute(rowQuery)
            print(rowQuery)
            cellDTSubset <- longTable@cells[eval(rowQuery), ]
        }
    #} else {
    #    cellDTSubset <- longTable@cells
    #}

    if (!missing(columnQuery)) {
        columnQuery <- substitute(columnQuery)
        print(columnQuery)
        drugDTSubset <- longTable@drugs[eval(columnQuery), ]
    } else {
        drugDTSubset <- longTable@drugs
    }

    rowIDs <- cellDTSubset$rowID
    colIDs <- drugDTSubset$colID

    if (missing(values)) { values <- names(longTable@dataList) }
    dataListIndexes <- values %in% names(longTable@dataList)

    valueData <- lapply(longTable@dataList[dataListIndexes],
                        FUN=.filterDataTable,
                        indexList=list(rowIDs, colIDs))

    return(LongTable(drugs=drugDTSubset, cells=cellDTSubset, dataList=valueData))
})

.filterDataTable <- function(DT, indexList) {
    if (length(indexList) > 2)
        stop("This object is 2D, please only pass in two ID vectors, one for
             rows and one for columns!")

    rowIndices <- indexList[[1]]
    colIndices <- indexList[[2]]

    return(copy(DT[rowID %in% rowIndices & colID %in% colIndices, ]))
}


#'
#'
#'
#'
#' @importFrom
setMethod('rowData')