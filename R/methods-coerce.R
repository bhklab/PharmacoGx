# ==== LongTable Class

#' Coerce a LongTable into a `data.table`
#'
#' @param from [`LongTable`] Object to coerce.
#' @param to [`character`] Class name to coerce to, currently only 'data.table'
#'    and 'data.frame' are supported
#'
#' return [`data.table`]
#'
#' @import data.table
#' @export
setAs('LongTable', 'data.table', def=function(from) {

    # local helpers
    .mapply <- function(...) mapply(..., SIMPLIFY=FALSE)

    # extract the assay data
    longTableData <- assays(from, key=TRUE)

    # join assays into a single table
    DT <- longTableData[[1]]
    longTableData[[1]] <- NULL
    for (i in seq_along(longTableData))
        DT <- merge.data.table(DT, longTableData[[i]], suffixes=c('', paste0('._', i)), by=.EACHI)

    # extract assay columns
    assayCols <- assayCols(from)

    # fix assayCols if there are duplicate column names between assays
    # the join will append '._n' where n is the assay index - 1
    ## TODO:: Make this a helper since it is reused in multiple functions
    .greplAny <- function(...) any(grepl(...))
    .paste0IfElse <- function(vector, suffix, isIn=c('rowKey', 'colKey'))
        ifelse(vector %in% isIn, vector, paste0(vector, suffix))
    hasSuffixes <- unlist(lapply(paste0('._', seq_along(longTableData)), .greplAny, x=colnames(DT)))
    if (any(hasSuffixes)) {
        whichHasSuffixes <- which(hasSuffixes) + 1
        assayCols[whichHasSuffixes] <-
            .mapply(FUN=.paste0IfElse,
                    vector=assayCols[whichHasSuffixes],
                    suffix=paste0('._', seq_along(longTableData))[hasSuffixes])
    }

    # join the row and column data
    DT <- merge.data.table(DT, rowData(longTable, key=TRUE))
    setkeyv(DT, c('rowKey', 'colKey'))
    DT <- merge.data.table(DT, colData(longTable, key=TRUE))

    # drop interal key columns
    DT[, c('rowKey', 'colKey') := NULL]

    # organize the returned columns
    colOrder <- c(setdiff(colnames(DT), unlist(assayCols)), unlist(assayCols))
    setcolorder(DT, colOrder)

    # capture configuration needed to reverse this operation
    ## TODO:: implement a configuration argument in constructor and get method for it.
    LongTable.config <- list(assayCols=assayCols,
                             rowDataCols=list(rowIDs(from), rowMeta(from)),
                             colDataCols=list(colIDs(from), colMeta(from)))
    attr(DT, 'LongTable.config') <- LongTable.config

    # return the data.table
    return(DT)
})

#' Coerce a LongTable into a `data.table`
#'
#' S3 version of coerce method for convenience.
#'
#' @param x [`LongTable`] to coerce to a `data.table`
#'
#' @return A [`data.table`] containing the data from the LongTable, as well
#'   as the `LongTable.config' attribute which contains the data needed to
#'   reverse the coercion.
#'
#' @export
as.data.table.LongTable <- function(x) as(x, 'data.table')

#' Coerce a LongTable into a `data.frame`
#'
#' Currently only supports coercing to data.table or data.frame
#'
#' @param from [`LongTable`] Object to coerce.
#' @param to [`character`] Class name to coerce to, currently only 'data.table'
#'    and 'data.frame' are supported
#'
#' @return [`data.table`] containing the data from the LongTable, with the
#'   `LongTable.config' attribute containg the metadata needed to reverse
#'    the coercing operation.
#'
#' @import data.table
#' @export
setAs('LongTable', 'data.frame', def=function(from){
    DT <- as(from, 'data.table')
    setDF(DT)
    return(DT)
})

#' Coerce a LongTable to a data.frame
#'
#' @param x [`LongTable`] to coerce to `data.frame`.
#'
#' @return [`data.frame`] containing the data from the LongTable, with the
#'   `LongTable.config' attribute containg the metadata needed to reverse
#'    the coercion operation.
#'
#' @export
as.data.frame.LongTable <- function(x) as(x, 'data.frame')


#' Coerce to data.table to LongTable
#'
#' Coerce a data.table with the proper configuration attributes back to a LongTable
#'
#' @param from A [`data.table`] with the 'LongTable.config' attribute, containing
#'    three lists named assayCols, rowDataCols and colDataCols. This attribute is
#'    automatically created when coercing from a LongTable to a data.table.
#'
#' @return [`LongTable`] object configured with the LongTable.config
#'
#' @export
setAs('data.table', 'LongTable', def=function(from) {

    if (!('LongTable.config' %in% names(attributes(from))))
        stop(.errorMsg('Coercing from data.table to LongTable only works if ',
                        'the LongTable.config attribute has been set.'))

    LongTable.config <- attr(from, 'LongTable.config')

    requiredConfig <- c('assayCols', 'rowDataCols', 'colDataCols')
    hasRequiredConfig <- requiredConfig %in% names(LongTable.config)
    if (!all(hasRequiredConfig))
        stop(.errorMsg('The LongTable.config attribute is missing the ',
            requiredConfig[!hasRequiedConfig], ' attribute(s).', collapse=', '))

    with(LongTable.config,
         buildLongTable(from, rowDataCols, colDataCols, assayCols))
})

#' Coerce to data.table to LongTable
#'
#' Coerce a data.table with the proper configuration attributes back to a LongTable
#'
#' @param x A [`data.table`] with the 'LongTable.config' attribute, containing
#'    three lists named assayCols, rowDataCols and colDataCols. This attribute is
#'    automatically created when coercing from a LongTable to a data.table.
#'
#' @return [`LongTable`] object configured with the LongTable.config
as.long.table <- function(x) as(x, 'LongTable')