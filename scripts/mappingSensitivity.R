# --- 0. Setup
library(data.table)

# define some tools for exploring dimensionality
length.unique <- function(...) length(unique(...))

setDTthreads(15)

# --- 0. Load Data

object <- .CCLE

raw <- as.data.table(sensitivityRaw(object))  # This automatically melts to long format
info <- as.data.table(sensitivityInfo(object), keep.rownames=TRUE)
profiles <- as.data.table(sensitivityProfiles(object), keep.rownames=TRUE)
info <- .info <- info[rn %in% unique(raw$V1)]
info[, experiment := seq_len(dim(.SD)[1]), by=.(cellid, drugid)]

# --- 1. Definite Mappings

mapped <- c('rn', 'drugid', 'cellid', 'experiment')

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
unmapped <- setdiff(colnames(info), mapped)


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

mapped <- c(mapped, possibleRow)

unmapped <- setdiff(colnames(info), mapped)
unmapped