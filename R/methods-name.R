#' name Getter
#'
#' Return the name of the PharmacoSet object
#'
#' @describeIn PharmacoSet Return the name of the PharmacoSet
#'
#' @examples
#' data(CCLEsmall)
#' name(CCLEsmall)
#'
#' @param object A [`PharmacoSet`]
#' @return `character` The name of the `PharmacoSet`
#'
#' @importMethodsFrom CoreGx name
#' @importFrom methods callNextMethod
#'
#' @export
setMethod('name', signature='PharmacoSet', function(object){
  callNextMethod(object=object)
})

#' name<- Setter
#'
#' Return the name of the `PharmacoSet` object
#'
#' @describeIn PharmacoSet Return the name of the PharmacoSet
#'
#'
#' @examples
#' data(CCLEsmall)
#' name(CCLEsmall) <- 'CCLEsmall'
#'
#' @param object A [`PharmacoSet`]
#' @return The name of the PharmacoSet
#'
#' @importFrom CoreGx name
#' @importFrom methods callNextMethod
#'
#' @importMethodsFrom CoreGx name<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod('name',
                 signature=signature(object='PharmacoSet', value='character'),
                 function(object, value){
  callNextMethod(object=object, value=value)
})
