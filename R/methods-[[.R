# ==== LongTable Class

#'
#'
#'
#' @describeIn LongTable
#'
#' @import crayon
#' @export
setMethod('[[', signature('LongTable'), function(x, i) {
    if (length(i) > 1)
        stop(cyan$bold("Please only select one assay! To subset on multiple assays please
            see ?subset"))
    else
        assay(x, i)
})