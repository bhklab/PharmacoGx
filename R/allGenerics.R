# ==== PharmacoSet accessor methods

# ---- drugInfo

#' drugInfo Generic
#'
#' Generic for drugInfo getter method
#'
#' @param object The \code{PharmacoSet} to retrieve drug info from
#'
#' @return A `data.frame` of annotations for drugs in the object
#'
#' @export
setGeneric("drugInfo", function(object) standardGeneric("drugInfo"))

#' drugInfo<- Generic
#'
#' Generic for drugInfo replace method
#'
#' @param object The [`PharmacoSet`] to replace drug info
#' @param value A `data.frame` with the new drug annotations
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
#' @return A `character` vector of drug names in the object
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

# ====