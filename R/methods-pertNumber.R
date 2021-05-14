#' pertNumber Getter
#'
#' Get the perturbation number for a specified `PharmcoSet` object
#'
#' @describeIn PharmacoSet Return the summary of available perturbation
#'   experiments
#'
#' @examples
#' data(CCLEsmall)
#' pertNumber(CCLEsmall)
#'
#' @param object A `PharmacoSet`
#'
#' @return A 3D  with the number of perturbation experiments per
#'   drug and cell line, and data type
#'
#' @importFrom CoreGx pertNumber
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(pertNumber, signature='PharmacoSet', function(object) {
    callNextMethod(object=object)
})


#' pertNumber<- Setter
#'
#' Set the perturbation number for a specified `PharmacoSet` object
#'
#' @examples
#' data(CCLEsmall)
#' pertNumber(CCLEsmall) <- pertNumber(CCLEsmall)
#'
#' @param object A \code{PharmacoSet}
#' @param value A new 3D \code{array} with the number of perturbation experiments
#'     per drug and cell line, and data type
#'
#' @return The updated `PharmacoSet`
#'
#' @describeIn PharmacoSet Update the summary of available perturbation
#'   experiments
#'
#' @importFrom CoreGx pertNumber<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod('pertNumber',
                 signature=signature(object='PharmacoSet', value='array'),
                 function(object, value) {
    callNextMethod(object=object, value=value)
})
