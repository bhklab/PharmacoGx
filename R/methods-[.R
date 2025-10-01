# ==== PharmacoSet Class
#'`[`
#'
#' @examples
#' data(CCLEsmall)
#' CCLEsmall["WM1799", "Sorafenib"]
#'
#' @param x object
#' @param i Cell lines to keep in object
#' @param j Drugs to keep in object
#' @param ... further arguments
#' @param drop A boolean flag of whether to drop single dimensions or not
#'
#'@return Returns the subsetted object
#'
#' @export
setMethod(`[`, 'PharmacoSet', function(x, i, j, ..., drop = FALSE) {
  validate_index <- function(idx, upper_bound, axis) {
    if (length(idx) == 0L) {
      return(integer(0))
    }
    if (anyNA(idx)) {
      stop(sprintf("Numeric %s indices must not contain NA values.", axis))
    }
    if (!is.numeric(idx)) {
      stop(sprintf("Numeric %s indices must be numeric.", axis))
    }
    if (any(idx <= 0)) {
      stop(sprintf("Numeric %s indices must be positive.", axis))
    }
    if (any(idx != round(idx))) {
      stop(sprintf("Numeric %s indices must be integers.", axis))
    }
    if (any(idx > upper_bound)) {
      stop(sprintf(
        "Numeric %s indices exceed dimension length (%d).",
        axis,
        upper_bound
      ))
    }
    as.integer(idx)
  }

  samples <- sampleNames(x)
  treatments <- treatmentNames(x)

  cells_requested_empty <- !missing(i) && length(i) == 0
  drugs_requested_empty <- !missing(j) && length(j) == 0

  if (missing(i)) {
    cell_names <- samples
  } else if (is.character(i)) {
    cell_names <- i
  } else if (is.numeric(i)) {
    idx <- validate_index(i, length(samples), "cell")
    cell_names <- if (length(idx)) samples[idx] else character(0)
    cells_requested_empty <- cells_requested_empty || (length(idx) == 0)
  } else {
    stop("Unsupported index type for cells; use numeric or character indices.")
  }
  if (!missing(i) && length(cell_names) == 0) {
    cells_requested_empty <- TRUE
  }

  if (missing(j)) {
    drug_names <- treatments
  } else if (is.character(j)) {
    drug_names <- j
  } else if (is.numeric(j)) {
    idx <- validate_index(j, length(treatments), "drug")
    drug_names <- if (length(idx)) treatments[idx] else character(0)
    drugs_requested_empty <- drugs_requested_empty || (length(idx) == 0)
  } else {
    stop("Unsupported index type for drugs; use numeric or character indices.")
  }
  if (!missing(j) && length(drug_names) == 0) {
    drugs_requested_empty <- TRUE
  }

  subsetTo(
    x,
    cells = cell_names,
    drugs = drug_names,
    molecular.data.cells = cell_names,
    allow.empty.cells = cells_requested_empty,
    allow.empty.drugs = drugs_requested_empty,
    ...
  )
})
