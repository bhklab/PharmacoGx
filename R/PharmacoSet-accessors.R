# Navigating this file:
# - Slot section names start with ----
# - Method section names start with ==
# 
# As a result, you can use Ctrl + f to find the slot or method you are looking
# for quickly, assuming you know its name.
# 
# For example Ctrl + f '== molecularProfiles' would take you the molecularProfiles
# method, while Ctrl +f '---- molecularProfiles' would take you to the slot
# section.


#' @include PharmacoSet-class.R allGenerics.R
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
## --------------

##
## == drugInfo

#' @rdname PharmacoSet-accessors
#'
#' @details
#' ## drug slot accessors
#'
#' __drugInfo__: `data.frame` Retrieve the drug metadata from a 
#' {`r .local_class`} objects `@drug` slot.
#'
#' @examples
#' drugInf <- drugInfo(CCLEsmall)
#'
#' @md
#' @aliases drugInfo{`r .local_class`}-method drugInfo
#' @exportMethod drugInfo
setMethod(drugInfo, signature="PharmacoSet", function(object) {
  object@drug
})

#' @rdname PharmacoSet-accessors
#'
#' @details
#' ## drug slot accessors
#'
#' __drugInfo__: `data.frame` Retrieve the drug metadata from a 
#' {`r .local_class`} objects `@drug` slot.
#'
#' @examples
#' drugInfo(`r .local_data`) <- drugInfo(`r .local_data`)
#'
#' @md
#' @aliases drugInfo<-{`r .local_class`},data.frame-method drugInfo<-
#' @exportMethod drugInfo<-
setReplaceMethod("drugInfo",
                 signature=signature(object="PharmacoSet", value="data.frame"),
                 function(object, value) {
  object@drug <- value
  object
})

## ==============
## ---- cell slot
## --------------

##
## == cellInfo

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_cellInfo(class_=.local_class)
#' @importFrom CoreGx cellInfo
#' @importFrom methods callNextMethod
setMethod(cellInfo, "PharmacoSet", function(object){
    callNextMethod(object)
})



#' @rdname PharmacoSet-accessors
#'
#' Generic for cellInfo replace method
#'
#' @examples
#' data(CCLEsmall)
#' cellInfo(CCLEsmall) <- cellInfo(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to replace cell info in
#' @param value A \code{data.frame} with the new cell annotations
#' @return Updated \code{PharmacoSet}
#'
#' @importFrom CoreGx cellInfo<-
#' @importFrom methods callNextMethod
setReplaceMethod("cellInfo", signature = signature(object="PharmacoSet", value="data.frame"), function(object, value){
  callNextMethod(object, value=value)
})