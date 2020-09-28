# ==== LongTable Class

#' Select an assay from a LongTable object
#'
#'
#'
#' @export
setMethod('$', signature('LongTable'), function(x, name) { x[[name]] })
