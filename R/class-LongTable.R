#' LongTable class for storing
#'
#'
#' @param drugs [`data.table`]
#' @param cells [`data.table`]
#' @param data [`data.table`]
#'
#' @return [`LongTable`] object containing the viability data from a treatment response experiment, indexes by a primary
#'  key
#'
#' @keywords internal
.LongTable <- setClass("LongTable",
                       slots=list(drugs="data.table",
                                  cells="data.table",
                                  data="data.table"))


#'
#'
#'
#'
#'
#'
LongTable <- function(drugs, cells, data) {

    if (!is(drugs, 'data.table')) {
        drugs <- data.table(drugs, keep.rownames='rownames')
    }

    if (!is(cells, 'data.table')) {
        cells <- data.table(cells, keep.rownames='rownames')
    }

    if (!is(data, 'data.table')) {
        data=data.table(data, keep.rownames=FALSE)
    }

    return(.LongTable(drugs=drugs, cells=cells, data=data))
}

dropCols <- 'combination_name'
cellCols <- c("cell_line", "BatchID")
drugCols <- c("drugA_name", "drugB_name")

#'
#'
#' @param filePath [`character`] Path to the .csv file containing the data and metadata from which to build the
#'  LongTable.
#' @param cellCols [`character`] Vector specifying one or more columns to be used for the cellData, corresponding to
#'  row subsetting for a `LongTable` object.
#' @param
#'
buildLongTableFromCSV <- function(filePath, cellCols, drugCols, dropCols) {

    rawData <- fread(filePath)

    drugData <- unique(rawData[, ..drugCols])[, colID := .I]
    cellData <- unique(rawData[, ..cellCols])[, rowID := .I]

    valueData <- rawData[drugData, on=drugCols][cellData, on=cellCols]
    valueData[, c(cellCols, drugCols, dropCols) := NULL]
    setkey(valueData, rowID, colID)
    setkey(drugData, colID)
    setkey(cellData, rowID)

    return(LongTable(drugs=drugData, cells=cellData, data=valueData))
}

setMethod('subset', signature('LongTable'), function())
