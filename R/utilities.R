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
.sensitivitySlotToLongTable <- function(object) {

    raw <- as.data.table(sensitivityRaw(object))  # This automatically melts to long format
                                                  # It also drops all NA rows
    info <- as.data.table(sensitivityInfo(object), keep.rownames=TRUE)
    profiles <- as.data.table(sensitivityProfiles(object), keep.rownames=TRUE)
    info[, drug_cell_rep := seq_len(dim(.SD)[1]), by=.(cellid, drugid)]

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
    sensDT <- merge.data.table(longRaw, info, all=TRUE)
    sensDT <- merge.data.table(sensDT, profiles, all=TRUE)

    # build identifiers assay identifiers
    assayRegex <- paste0(assayIDs, '_\\d+')
    .getAssayCols <- function(name, colNames) grep(name, colNames, value=TRUE)
    assayCols <- lapply(assayRegex, .getAssayCols, colNames=colnames(sensDT))
    names(assayCols) <- tolower(assayIDs)

    mappings <- .getLongTableDimensionMappings(info)

    rowDataCols <- list(c('cellid', 'drug_cell_rep'), mappings$rowMeta)
    colDataCols <- list(c('drugid', 'drug_cell_rep'), mappings$colMeta)

    if (length(mappings$unmapped) > 1)
        assayCols$experiment_metadata <- mappings$unmapped

    assayCols$sensitivity_profiles <- colnames(profiles[, -'rn'])

    ## FIXME:: Re-implement metadata once we sort out column mappings.
    #metadata <- as.list(unique(sensDT[, mappings$metadata, with=FALSE]))

    longTable <- buildLongTable(from=sensDT, rowDataCols, colDataCols, assayCols)

    return(longTable)
}

## TODO:: refactor this to be shorter/more concise
#'
#'
#' @import data.table
#' @export
.getLongTableDimensionMappings <- function(info) {

    # define some tools for exploring dimensionality
    length.unique <- function(...) length(unique(...))

    mapped <- c('rn', 'drugid', 'cellid', 'drug_cell_rep')

    # Determine which columns map to metadata
    # (i.e., which columns are the same)
    byMeta <- sapply(info, length.unique)
    mapsToMeta <- byMeta == 1
    definiteMeta <- names(which(mapsToMeta))

    mapped <- c(mapped, definiteMeta)

    # Determine which columns map to rn
    # (i.e., which columns are all unique)
    byRow <- sapply(info, length.unique)
    mapsToRow <- byRow == nrow(info)
    definiteRow <- names(mapsToRow[mapsToRow])
    definiteRow <- setdiff(definiteRow, mapped)

    mapped <- c(mapped, definiteRow)

    # Determine which columns map to drug
    byDrug <- info[, lapply(.SD, length.unique), by=drugid]
    maxByDrug <- sapply(byDrug[, -'drugid'], max)
    definiteDrug <- names(which(maxByDrug == 1))
    definiteDrug <- setdiff(definiteDrug, mapped)

    mapped <- c(mapped, definiteDrug)

    # Determine which columns map to cell
    byCell <- info[, lapply(.SD, length.unique), by=cellid]
    maxByCell <- sapply(byCell[, -'cellid'], max)
    definiteCell <- names(which(maxByCell == 1))
    definiteCell <- setdiff(definiteCell, mapped)

    mapped <- c(mapped, definiteCell)

    # ---- 2. Possible Mappings

    # Possible metadata
    maybeMeta <- names(which(byMeta < 5))
    maybeMeta <- setdiff(maybeMeta, mapped)

    mapped <- c(mapped, maybeMeta)

    # Possible cell mappings
    maybeCell <- names(which(maxByCell < 5))
    maybeCell <- setdiff(maybeCell, mapped)

    # Possible drug mappings
    maybeDrug <- names(which(maxByDrug < 5))
    maybeDrug <- setdiff(maybeDrug, mapped)

    # Overlap in possible cell and drug mappings
    maybeBoth <- setdiff(intersect(maybeCell, maybeDrug), maybeMeta)
    if (length(maybeBoth) > 0)
        stop("Need to handle ambiguous mapping between cell/drug")
        betterDrug <- maxByDrug[maybeBoth] <= maxByCell[maybeBoth]
        ## TODO:: Handle this if a case occurs

    definiteMeta <- c(definiteMeta, maybeMeta)
    definiteDrug <- c(definiteDrug, maybeDrug)
    definiteCell <- c(definiteCell, maybeCell)

    mapped <- c(mapped, maybeDrug, maybeCell)

    # Possible row mappings
    maxDrugCell <- max(byRow[c('cellid', 'drugid')])
    possibleRow <- names(which(byRow < maxDrugCell & byRow != nrow(info)))
    possibleRow <- setdiff(possibleRow, mapped)

    definiteRow <- c(definiteRow, possibleRow)

    mapped <- c(mapped, possibleRow)

    unmapped <- setdiff(colnames(info), mapped)
    unmapped <- unique(c('rn', definiteRow, definiteMeta, unmapped))

    ## FIXME:: Removed metadata columns to ensure I can recreate sensitivityInfo
    mappings <- list(rowMeta=definiteCell, colMeta=definiteDrug, unmapped=unmapped)
    return(mappings)
}