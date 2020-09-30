# ==== PharmacoSet accessor methods

# ---- drugInfo

#' drugInfo Generic
#'
#' Generic for drugInfo getter method
#'
#' @param object The \code{PharmacoSet} to retrieve drug info from
#'
#' @return A [`data.frame`] of annotations for drugs in the object
#'
#' @export
setGeneric("drugInfo", function(object) standardGeneric("drugInfo"))

#' drugInfo<- Generic
#'
#' Generic for drugInfo replace method
#'
#' @param object The [`PharmacoSet`] to replace drug info
#' @param value A [`data.frame`] with the new drug annotations
#'
#' @return The [`object`] with updated drug annotations
#'
#' @export
setGeneric('drugInfo<-', function(object, value) standardGeneric('drugInfo<-'))

# --- drugNames

#' drugNames Generic
#'
#' A generic for the drugNames method
#'
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall)
#'
#' @param object The [`PharmacoSet`] to return drug names from
#'
#' @return A [`character`] vector of drug names in the object
#'
#' @export
setGeneric("drugNames", function(object) standardGeneric("drugNames"))

#' drugNames<- Generic
#'
#' A generic for the drugNames replacement method
#'
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall) <- drugNames(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{character} vector of the new drug names
#'
#' @return The [`object`] with updated drug names
#'
#' @export
setGeneric("drugNames<-", function(object, value) standardGeneric("drugNames<-"))

# ==== LongTable Class

#' Generic method for resetting indexing in an S4 object
#'
#' This method allows integer indexes used to maintain referential integrity
#'   internal to an S4 object to be reset. This is useful particularly after
#'   subsetting an object, as certain indexes may no longer be present in the
#'   object data. Reindexing removes gaps integer indexes and ensures that the
#'   smallest contiguous integer values are used in an objects indexes.
#'
#' @param object [`S4`] An object to redo indexing for
#' @param ... [`pairlist`] Allow definition of new parameters to this generic.
#'
#' @export
setGeneric('reindex', function(object, ...) standardGeneric('reindex'))

#' Build a LongTable object
#'
#' @param from What to build the LongTable from?
#' @param ... [`pairlist`] Allow definition of new parameters for
#'     implementations of this generic.
#'
#' @export
setGeneric('buildLongTable', function(from, ...) standardGeneric('buildLongTable'))


# ===== Other Generics

# ---- plot

# FIXME:: Surely we can import this generic from somewhere?
setGeneric("plot")


