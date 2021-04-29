#' @include class-PharmacoSet.R allGenerics.R
NULL

.local_class <- 'PharmacoSet'

# =======================================
# Accessor Method Documentation Object
# ---------------------------------------

#' @eval CoreGx:::.docs_CoreSet_accessors(name_='PharmacoSet-accessors', 
#'    class_='PharmacoSet')
NULL


# ======================================
# Accessor Methods
# --------------------------------------


## --------------
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