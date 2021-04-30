#' @include class-PharmacoSet.R allGenerics.R
NULL

.local_class <- 'PharmacoSet'
.local_data <- 'CCLEsmall'
data(CCLEsmall)

# =======================================
# Accessor Method Documentation Object
# ---------------------------------------

#' @name PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_accessors(class_=.local_class)
NULL


# ======================================
# Accessor Methods
# --------------------------------------



## ==============
## ---- drug slot

#' @examples
#' data(CCLEsmall)
#' drugInf <- drugInfo(CCLEsmall)
#'
#' @rdname PharmacoSet-accessors
#'
#' @md
#' @export
setMethod(drugInfo, signature="PharmacoSet", function(object) {
  object@drug
})

#' @examples
#' data(CCLEsmall)
#' drugInf <- drugInfo(CCLEsmall)
#'
#' @rdname PharmacoSet-accessors
#'
#' @md
#' @export
setReplaceMethod("drugInfo",
                 signature=signature(object="PharmacoSet", value="data.frame"),
                 function(object, value) {
  object@drug <- value
  object
})

## ==============
## ---- cell slot

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_cellInfo(class_=.local_class)
#' @importFrom CoreGx cellInfo
#' @importFrom methods callNextMethod
setMethod(cellInfo, "PharmacoSet", function(object){
  callNextMethod(object)
})

