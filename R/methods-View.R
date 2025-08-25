#' Guarded View for CoreSet/PharmacoSet
#'
#' Emits the same compatibility error as `CoreGx::show()` and stops for
#' outdated CoreSet-derived objects. Otherwise delegates to whatever
#' `utils::View` is installed in `package:utils` (RStudio/Positron override
#' or base).
#'
#' @param x Any R object, typically a data.frame-like or `CoreSet`-derived.
#' @param title Optional title for the data viewer.
#'
#' @return Invisibly returns whatever the underlying viewer returns.
#' @seealso [utils::View()]
#'
#' @examples
#' if (interactive()) {
#'   data("CCLEsmall", package = "PharmacoGx")
#'   View(CCLEsmall)
#' }
#'
#' @export
View <- function(x, title = NULL) {
  # Only intercept for CoreSet-derived objects; keep everything else untouched
  if (inherits(x, "CoreSet")) {
    # Replicate CoreGx::show() outdated check (CoreGx/R/CoreSet-class.R:424-428)
    sn <- methods::slotNames(x)
    hasSample <- "sample" %in% sn
    hasTreatment <- "treatment" %in% sn
    if (!(hasSample && hasTreatment)) {
      # Hard stop with the same message as CoreGx::show()
      stop(
        CoreGx::.errorMsg(
          "This ",
          class(x)[1],
          " object appears to be out of date! ",
          "Please run object <- updateObject(object) to update ",
          "the object for compatibility with the current release."
        ),
        call. = FALSE
      )
    }
  }

  # Ensure Positron/RStudio/base receive a scalar string title
  if (missing(title) || is.null(title)) {
    title2 <- base::deparse(substitute(x))
  } else {
    title2 <- base::as.character(title)[1]
    if (is.na(title2)) title2 <- base::deparse(substitute(x))
  }

  # Delegate to whatever View is installed in package:utils
  get("View", envir = as.environment("package:utils"))(x, title2)
}
