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

# ==== LongTable Accessor Methods

# ---- Private Helper Methods

#' Return the identifiers for the column meta data in an object
#'
#' @export
#' @keywords internal
setGeneric('.colIDData', function(object, ...) standardGeneric('.colIDData'))

#' Return the identifiers for the row metadata columns in an object
#'
#'
#' @export
#' @keywords
setGeneric('.rowIDData', function(object, ...) standardGeneric('.rowIDData'))


#' Private method to retrieve the .colIDs property from an object
#'
#' @export
#' @keywords internal
setGeneric('.colIDs', function(object, ...) standardGeneric('.colIDs'))

#' Private method to retrieve the .rowIDs property from an object
#'
#' @export
#' @keywords internal
setGeneric('.rowIDs', function(object, ...) standardGeneric('.rowIDs'))

#' Private method to retrieve the both the .rowIDs and .colIDs properties from
#'   an object in a list.
#'
#' @param object [`any`] An object with a .intern slot with the items .rowIDs
#'   and .colIDs
#'
#' @return A [`list`] v
#'
#' @export
#' @keywords internal
setGeneric('.dimIDs', function(object, ...) standardGeneric('.dimIDs'))



# ===== Other Generics

# ---- plot

# FIXME:: Surely we can import this generic from somewhere?
setGeneric("plot")


