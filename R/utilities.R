#' Convenience function for collapsing a character vector
#'
#' @param ... [`pairlist`] One or more character vectors
#' @param collapse [`character`] Argument to collapse of paste0, default is ' '.
#'
#' @return [`character`] A single character vector.
#'
#' @keywords internal
#' @export
#' @noRd
.collapse <- function(..., collapse=' ')
    paste0(..., collapse=collapse)

#' Returns a colorized error message (magenta)
#'
#' @param ... [`pairlist`] One or more strings or character vectors, also
#'   accepts any params to paste0.
#'
#' @return [`character`] Colorized string with results from paste0(...)
#'
#' @keywords internal
#' @export
#' @noRd
.errorMsg <- function(...) magenta$bold(paste0(...))

#' Returns a colorized warning message (cyan)
#'
#' @param ... [`pairlist`] One or more strings or character vectors, also
#'   accepts any params to paste0.
#'
#' @return [`character`] Colorized string with results from paste0(...)
#'
#' @keywords internal
#' @export
#' @noRd
.warnMsg <- function(...) cyan$bold(paste0(...))

#' Get the types of all items in a list
#'
#' @param list A [`list`] to get the types from
#' @param ... [`pairlist`] Additional arguments to FUN
#' @param FUN [`function`] or [`character`] Either a function, or the name
#'   of a function which returns a single logical value. The default function
#'   uses `is`, specify the desired type in `...`. You can also use other
#'   type checking functions such as is.character, is.numeric, or is.data.frame.
#'
#' @export
#' @noRd
is.items <- function(list, ..., FUN=is)
    vapply(list, FUN=FUN, FUN.VALUE=logical(1), ...)

#'
#'
#'
#' @import data.table
#' @export
sensitivitySlotToLongTable <- function(object) {

    raw <- as.data.table(sensitivityRaw(object))  # This automatically melts to long format
    info <- as.data.table(sensitivityInfo(object), keep.rownames=TRUE)
    profiles <- as.data.table(sensitivityProfiles(object), keep.rownames=TRUE)

    # preprocess raw array
    setnames(raw, seq_len(3), c('rn', 'replicate', 'assay'))
    assayIDs <- unique(raw$assay)
    raw[, value := as.numeric(value)]
    raw[, replicate := as.integer(gsub('\\D*', '', replicate))]
    # Split value into one column for each assay (long -> wide)
    longRaw <- dcast(raw, rn + replicate ~ ..., value.var=c('value'))
    # Split assay columns into assay per replicate (long -> wide)
    longRaw <- dcast(longRaw, rn ~ replicate, value.var=assayIDs)

    # set keys for joins
    setkeyv(info, 'rn')
    setkeyv(profiles, 'rn')
    setkeyv(longRaw, 'rn')

    # join raw and info to the assay data
    sensDT <- merge.data.table(longRaw, info)
    sensDT <- merge.data.table(sensDT, profiles)

    ## TODO:: Implement dimensionality check to ensure no rows are being dropped in joins

    # build identifiers assay identifiers
    .getAssayCols <- function(name, colNames) grep(name, colNames, value=TRUE)
    assayCols <- lapply(assayIDs, .getAssayCols, colNames=colnames(sensDT))
    names(assayCols) <- tolower(assayIDs)

    # find columns which match 1:1 with rownames
    .length.unique <- function(...) length(unique(...))
    uniquelyMapsToRn <- vapply(info, .length.unique, integer(1)) == nrow(info)
    infoNoRn <- info[, .SD, .SDcols=!uniquelyMapsToRn]

    # build row identifiers
    .length.unique.is.1 <- function(...) length(unique(...)) == 1
    uniquelyMapsToCellid <- vapply(infoNoRn[, lapply(.SD, .length.unique.is.1), by=cellid][, -c('cellid')], all, logical(1))
    rowMeta <- names(uniquelyMapsToCellid[uniquelyMapsToCellid])

    # build column identifiers
    uniquelyMapsToDrugid <- vapply(infoNoRn[, lapply(.SD, .length.unique.is.1), by=drugid][, -c('drugid')], all, logical(1))
    colMeta <- names(uniquelyMapsToDrugid[uniquelyMapsToDrugid])

    # get columns for metadata
    metadataCols <- intersect(rowMeta, colMeta)
    metadata <- as.list(unique(infoNoRn[, ..metadataCols]))
    metadata <- lapply(metadata, unique)

    # remove metadata colums from row and column meta
    rowMeta <- setdiff(rowMeta, metadataCols)
    colMeta <- setdiff(colMeta, metadataCols)

    # find the minimum columns in addition to the id columns to uniquely identify each row
    ## TODO:: Refactor this to something not gross looking
    if (nrow(unique(info[, .(drugid, cellid)])) != nrow(info)) {
        candidateCols <- setdiff(colnames(infoNoRn), unique(c(rowMeta, colMeta, metadataCols, 'drugid', 'cellid')))
        uniquelyIdentify <- FALSE
        i <- 0
        idCols <- list()
        while (!uniquelyIdentify) {
            i <- i + 1
            print(i)
            colCombos <- if (i == 1) candidateCols else combn(candidateCols, i, simplify=FALSE)
            colCombos <- lapply(colCombos, c, c('drugid', 'cellid'))
            for (columns in colCombos) {
                if (nrow(info) == nrow(unique(info[, ..columns]))) {
                    print(idCols)
                    idCols <- c(idCols, list(columns))
                    uniquelyIdentify <- TRUE
                } else {
                    if (i == length(candidateCols)) {
                        idCols <- NULL
                        uniquelyIdentify <- TRUE
                    }
                }
            }
        }
    }

    # error if can't find id columns
    if (is.null(idCols))
        stop(.errorMsg('Unable to find columns uniquely identifying each row.',
             'Please try manuallly building the LongTable.'))

    # if there are more than one matches, determine which match results in less unique rows total
    # (this should help minimize the dimensionality of the LongTable)
    additionalIDCols <- lapply(idCols, setdiff, c('drugid', 'cellid'))
    sumUniqueCols <- c()
    for (columns in additionalIDCols) {
        sumUniqueCols <- c(sumUniqueCols,
            sum(info[, lapply(.SD, .length.unique), .SDcols=columns]))
    }
    idCols <- idCols[[which.min(sumUniqueCols)]]
    additionalIDCols <- additionalIDCols[[which.min(sumUniqueCols)]]

    # determine whether the additional id columns should go into sample or drug info
    additionalColumnCombos <- lapply(seq_along(additionalIDCols), combn, x=additionalIDCols, simplify=FALSE)
    if (is.list(additionalColumnCombos[[1]]))
        additionalColumnCombos <- unlist(additionalColumnCombos, recursive=FALSE)
    additionalColumnCombos <- c(additionalColumnCombos, list(character(0)))
    otherColumns <- lapply(additionalColumnCombos, setdiff, x=additionalIDCols)

    longTableDims <- list()  # Note: Not really dims, was getting integer overflow when multipling...
    for (i in seq_along(otherColumns)) {
        longTableDims[[i]] <- sum(
            nrow(unique(info[, .SD, .SDcols=c('cellid', otherColumns[[i]])])),
            nrow(unique(info[, .SD, .SDcols=c('drugid', additionalColumnCombos[[i]])]))
        )
    }
    rowIDs <- c('cellid', additionalColumnCombos[[which.min(longTableDims)]])
    colIDs <- c('drugid', otherColumns[[which.min(longTableDims)]])

    # build the LongTable object
    unmappedCols <- setdiff(candidateCols, unique(unlist(c(rowIDs, rowMeta, colIDs, colMeta, metadataCols))))

    # put unmapped columns in an assay
    assayCols <- c(assayCols, list(unmapped=unmappedCols))

    # assemble row and column data
    rowDataCols <- list(rowIDs, rowMeta)
    colDataCols <- list(colIDs, colMeta)

    longTable <- buildLongTable(from=sensDT, rowDataCols, colDataCols, assayCols)
    ## TODO:: Add metadata parameter to buildLongTable function
    metadata(longTable) <- metadata

    return(longTable)
}