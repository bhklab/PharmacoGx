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

    # validate input and return useful messages if invalid
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
    rowDataColnames <- lapply(rowDataCols, names)
    notNullRownames <- !vapply(rowDataColnames, FUN=is.null, FUN.VALUE=logical(1))
    if (any(notNullRownames))
        for (i in which(notNullRownames)) {
            setnames(rowData, rowDataCols[[i]], names(rowDataCols[[i]]))
            rowDataCols[[i]] <- names(rowDataCols[[i]])
        }

    colDataColnames <- lapply(colDataCols, names)
    notNullColnames <- !vapply(colDataColnames, FUN=is.null, FUN.VALUE=logical(1))
    if (any(notNullColnames))
        for (i in which(notNullColnames)) {
            setnames(colData, colDataCols[[i]], names(colDataCols[[i]]))
            colDataCols[[i]] <- names(colDataCols[[i]])
        }

    # add the index columns to the different assay column vectors
    # this allows the .selectDataTable helper to be more general
    .prependToVector <- function(vector, values) c(values, vector)
    assayCols <- lapply(assayCols, FUN=.prependToVector, values=c('rowKey', 'colKey'))
    if (is.null(names(assayCols))) names(assayCols) <- paste0('assay', seq_along(assayCols))
    assays <- lapply(assayCols, .selectDataTable, DT=assayData)

    # remove the colname suffixes by reference from assays which had the same
    # colnames prior to joining into a single DT
    for (assay in assays) {
        setnames(assay, colnames(assay), gsub('\\._\\d$', '', colnames(assay)))
    }

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
    if (length(from) > 1)  # Call list subsetting method
        buildLongTable(as.list(from), rowDataCols, colDataCols, assayCols)
    if (!file.exists(from))
        stop(magenta$bold("The is no file at path: ", from, '. Please double
            check the location of the source file!'))

    # read in data
    tableData <- .freadNA(filePath)

    return(buildLongTable(from=tableData, rowDataCols, colDataCols, assayCols))
})

#' Create a LongTable object from a list containing file paths, data.frames and
#'   data.tables.
#'
#' @param from [`list`] A list containing any combination of character file paths,
#'  data.tables and data.frames which will be used to construct the LongTable.
#'
#'
#' @import data.table
#' @importFrom crayon magenta cyan
#' @export
setMethod('buildLongTable', signature(from='list'),
          function(from, rowDataCols, colDataCols, assayCols) {

    # local helpers
    .mapply <- function(...) mapply(..., SIMPLIFY=FALSE)

    # preprocess from list
    isChar <- is.items(from, 'character')
    isDT <- is.items(from, FUN=is.data.table)
    isDF <- is.items(from, FUN=is.data.frame) & !isDT

    if (!all(isChar | isDT | isDF))
        stop(.errorMsg('List items at indexes ',
             which(!(isChar | isDT | isDF )),
             ' are not character, data.table or data.frame.', collapse=', '))

    if (any(isChar)) from <- c(from[!isChar], lapply(from[isChar], FUN=.freadNA))
    if (any(isDF)) for (i in which(isDF)) setDT(from[[i]])

    # validate mappings
    idCols <- c(unlist(rowDataCols[[1]]), unlist(colDataCols[[1]]))
    dataColNames <- lapply(from, FUN=colnames)
    idColsIn <- lapply(dataColNames, `%in%`, x=idCols)
    hasAllIdCols <- unlist(lapply(idColsIn, FUN=all))
    if (!all(hasAllIdCols)) {
        missingCols <- unique(unlist(.mapply(`[`, x=idCols, i=idColsIn)))
        stop(.errorMsg('Assays ', which(hasAllIdCols),
             ' are missing one or more id columns: ', missingCols,
             collapse=', '))
    }

    # join assays into a single table
    DT <- from[[1]]
    from[[1]] <- NULL
    assayNames <- names(assayCols)
    for (i in seq_along(from))
        DT <- merge.data.table(DT, from[[i]], on=idCols, all=TRUE,
            suffixes=c('', paste0('._', i)))

    # fix assayCols if there are duplicate column names between assays
    # the join will append '._n' where n is the assay index - 1
    .greplAny <- function(...) any(grepl(...))
    hasSuffixes <- unlist(lapply(paste0('._', seq_along(from)), .greplAny, x=colnames(DT)))
    if (any(hasSuffixes)) {
        whichHasSuffixes <- which(hasSuffixes) + 1
        assayCols[whichHasSuffixes] <-
            .mapply(paste0, assayCols[whichHasSuffixes],
                    paste0('._', seq_along(from))[hasSuffixes])
    }

    # construct new LongTable
    buildLongTable(from=DT, rowDataCols, colDataCols, assayCols)

})


# ---- Helper Methods

#' fread with more default na.strings
#'
#' @keywords internal
#' @export
#' @noRd
.freadNA <- function(...) {
    as.na <- unique(c(getOption('datatable.na.string'),
        c('NA', 'NULL', 'NaN', 'missing', 'None',
            'none', 'na', 'null', 'Null', 'Na')))
    fread(..., na.strings=as.na)
}


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
#' @export
#' @noRd
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