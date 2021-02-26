#' sensNumber Getter
#'
#' Get the sensitivity numbers for a `PharmacoSet` object
#'
#' @describeIn PharmacoSet Return the summary of available sensitivity
#'   experiments
#'
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall)
#'
#' @param object A \code{PharmacoSet}
#' @return A \code{data.frame} with the number of sensitivity experiments per
#'   drug and cell line
#'

#'
#' @importFrom CoreGx sensNumber
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(sensNumber, 'PharmacoSet', function(object){
  callNextMethod(object=object)
})

#' Return the number of times each row-column combination occurs in a LongTable
#'
#' Reconstruct the @sensitivity$n list data from a LongTable object. This allows
#'     backwards compatibility with the current accessors for the @sensitivity
#'     list object.
#'
#' @section WARNING:
#' Because a LongTable has incomplete information about the rows
#'    and columns present in a CoreSet, this function is unable to zero
#'    pad missing rows and columns. This will need to be implemented in the
#'    sensNumber method for a class inheriting from the CoreSet. For example,
#'   in a `PharmacoGx::PharmacoSet`, `LongTable` rows are cells and columns
#'   are drugs.
#'
#' @param longTable A [`LongTable`] longTable object to rebuild the results
#'   of sensNumber for.
#'
#' @return [`matrix`] A row by column matrix containing a count for the number
#'   of times a row-column combination occurs in a LongTable object.
#'
#' @keywords internal
#' @noRd
.rebuildN <- function(longTable) {

    # Extract the information needed to reconstruct the sensitivityRaw array
    meta <- assay(longTable, 'experiment_metadata')[, .(rn, rowKey, colKey)]
    setkeyv(meta, c('rowKey', 'colKey'))
    rowData <- rowData(longTable, key=TRUE)[, .(cellid, drug_cell_rep, rowKey)]
    setkeyv(rowData, 'rowKey')
    colData <- colData(longTable, key=TRUE)[, .(drugid, drug_cell_rep, colKey)]
    setkeyv(colData, 'colKey')

    # join the tables into the original data
    num <- merge.data.table(meta, rowData, all=TRUE)
    setkeyv(num, 'colKey')
    num <- merge.data.table(num, colData, all=TRUE)[, .(cellid, drugid, drug_cell_rep.x)]

    num <- dcast(num, cellid ~ drugid, value.var='drug_cell_rep.x',
                 fun.aggregate=max, fill=0)
    #num <- unique(num)
    rownames <- num$cellid
    num[, cellid := NULL]

    num <- as.matrix(num)
    rownames(num) <- rownames
    return(num)
}

## TODO:: Make this a unit test
## all.equal(num[rownames(SN), colnames(SN)], SN)