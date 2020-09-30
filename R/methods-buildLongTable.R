# ==== LongTable Class

#' Create a LongTable object from a single data.table or data.frame object
#'
#' @param from [`character`] Path to the .csv file containing the data and
#'   metadata from which to build the `LongTable`.
#' @param colDataCols [`list`] List with two `character` vectors, the first
#'   specifying one or more columns to be used as column identifiers (e.g.,
#'   drug name columns) and the second containing any additional metadata
#'   columns related to the column identifiers. If you wish to rename any of
#'   these columns, assign the new names to their respective character vectors.
#' @param rowDataCols [`list`] List with two `character` vectors, the first
#'   specifying one or more columns to be used as cell identifiers (e.g.,
#'   cell-line name columns) and the second containing any additional metadata
#'   columns related to the cell identifiers. If you wish to rename any of
#'   these columns, assign the new names to their respective character vectors.
#' @param assayCols [`list`] A named list of character vectors specifying how to
#'   parse assay columns into a list of `data.table`s. Each list data.table
#'   will be named for the name of corresponding list item and contain the columns
#'   specified in the character vector of column names in each list item. If
#'   there are no names for assayCols, the assays will be numbered by instead.
#'
#' @return A [`LongTable`] object containing one or more assays, indexed by
#'   rowID and colID.
#'
#' @import data.table
#' @export
setMethod('buildLongTable', signature(from='data.frame'),
          function(from, rowDataCols, colDataCols, assayCols)
{
    # handle missing params
    missingParams <- c(missing(rowDataCols), missing(colDataCols), missing(assayCols))
    if (any(missingParams))
        stop(magenta$bold('The following parameters are required:',
            c('rowDataCols', 'colDataCols', 'assayCols')[missingParams]))

    ## TODO:: Check input parameters are valid

    # convert to data.table
    if (!is.data.table(from))
        from <- data.table(from, keep.rownames=FALSE)

    # build drug and cell metadata tables and index by the appropriate ID
    colData <- unique(from[, unlist(colDataCols), with=FALSE])
    colData[, colKey := seq_len(.N)]
    rowData <- unique(from[, unlist(rowDataCols), with=FALSE])
    rowData[, rowKey := seq_len(.N)]

    # add the row and column ids to the value data; as.character strips names
    # which were causing an error inside data.table call
    assayData <- from[rowData, on=as.character(unlist(rowDataCols))][colData, on=as.character(unlist(colDataCols))]
    rm(from)
    assayData[, as.character(c(unlist(rowDataCols), unlist(colDataCols))) := NULL]
    setkey(assayData, rowKey, colKey)

    setkey(rowData, rowKey)
    setkey(colData, colKey)

    # rename columns, if necessary
    ## TODO:: Refactor this mess
    rowDataColnames <- lapply(rowDataCols, names)
    if (length(rowDataColnames) < 2) rowDataColnames <- list(rowDataColnames)
    notNullRownames <- !vapply(rowDataColnames, FUN=is.null, FUN.VALUE=logical(1))
    if (any(notNullRownames)) {
        for (i in which(notNullRownames)) {
            setnames(rowData, rowDataCols[[i]], names(rowDataCols[[i]]))
        }
        rowDataCols[[i]] <- names(rowDataCols[[i]])
    }
    colDataColnames <- lapply(colDataCols, names)
    if (length(colDataColnames) < 2) colDataColnames <- list(colDataColnames)
    notNullColnames <- !vapply(colDataColnames, FUN=is.null, FUN.VALUE=logical(1))
    if (any(notNullColnames)) {
        for (i in which(notNullColnames)) {
            setnames(colData, colDataCols[[i]], names(colDataCols[[i]]))
        }
        colDataCols[[i]] <- names(colDataCols[[i]])
    }

    # add the index columns to the different assay column vectors
    # this allows the .selectDataTable helper to be more general
    .prependToVector <- function(vector, values) c(values, vector)
    assayCols <- lapply(assayCols, FUN=.prependToVector, values=c('rowKey', 'colKey'))
    if (is.null(names(assayCols))) names(assayCols) <- paste0('assay', seq_along(assayCols))
    assays <- lapply(assayCols, .selectDataTable, DT=assayData)

    ## applicable to any type of data. Maybe allow user to specify? For example
    ## by naming the elements of rowDataCols and colDataCols?
    return(LongTable(rowData=rowData, rowIDs=rowDataCols[[1]],
                     colData=colData, colIDs=colDataCols[[1]],
                     assays=assays))
})

#' Create a LongTable object from a single .csv file
#'
#' @param from [`character`] Path to the .csv file containing the data and
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
#' @importFrom crayon magenta
#' @export
setMethod('buildLongTable', signature(from='character'),
          function(from, rowDataCols, colDataCols, assayCols)
{
    if (!file.exists(from))
        stop(magenta$bold("The is no file at path: ", from, '. Please double
            check the location of the source file!'))

    # read in data
    tableData <- fread(filePath,
                na.strings=unique(c(getOption('datatable.na.string'),
                                  c('NA', 'NULL', 'NaN', 'missing', 'None',
                                    'none', 'na', 'null', 'Null', 'Na'))))

    return(buildLongTable(from=tableData, rowDataCols, colDataCols, assayCols))
})

## TODO:: Implement buildLongTable method for lists of objects or lists of file paths.



# ---- Helper Methods

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