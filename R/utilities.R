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