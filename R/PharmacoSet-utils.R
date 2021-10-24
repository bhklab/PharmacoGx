#' @include PharmacoSet-class.R PharmacoSet-accessors.R
NULL

.local_class <- 'PharmacoSet'
.local_data <- 'CCLEsmall'
.local_treatment <- 'drug'


#### PharmacoGx dynamic documentation
####
#### Warning: for dynamic docs to work, you must set
#### Roxygen: list(markdown=TRUE, r6=FALSE)
#### in the DESCRPTION file!


# ===================================
# Utility Method Documentation Object
# -----------------------------------


#' @name PharmacoSet-utils
#' @eval CoreGx:::.docs_CoreSet_utils(class_=.local_class)
#' @eval .parseToRoxygen("@examples data({data_})", data_=.local_data)
NULL


# ======================================
# Subset Methods
# --------------------------------------



## ===================
## ---- subsetBySample
## -------------------


#' @rdname PharmacoSet-utils
#' @importMethodsFrom CoreGx subsetBySample
#' @eval CoreGx:::.docs_CoreSet_subsetBySample(class_=.local_class,
#' data_=.local_data)
setMethod('subsetBySample', signature(x='PharmacoSet'), function(x, samples) {
    callNextMethod(x=x, samples=samples)
})


## ======================
## ---- subsetByTreatment
## ----------------------


#' @rdname PharmacoSet-utils
#' @importMethodsFrom CoreGx subsetByTreatment
#' @eval CoreGx:::.docs_CoreSet_subsetByTreatment(class_=.local_class, 
#' data_=.local_data, treatment_=.local_treatment)
setMethod('subsetByTreatment', signature(x='PharmacoSet'),
        function(x, treatments) {
    callNextMethod(x=x, treatments=treatments)
})


## ====================
## ---- subsetByFeature
## --------------------


#' @rdname PharmacoSet-utils
#' @importFrom CoreGx subsetByFeature
#' @eval CoreGx:::.docs_CoreSet_subsetByFeature(class_=.local_class, 
#' data_=.local_data)
setMethod('subsetByFeature', signature(x='PharmacoSet'), 
        function(x, features, mDataTypes) {
    callNextMethod(x=x, features=features, mDataTypes)
})

## ===========
## ---- subset
## -----------

#'
#' 
#' 
setMethod('subset', signature('PharmacoSet'),
        function(x, samples, treatments, features, ..., mDataTypes) {
    callNextMethod(x=x, samples=samples, treatments=treatments, 
        features=features, ..., mDataTypes=mDataTypes)
})
