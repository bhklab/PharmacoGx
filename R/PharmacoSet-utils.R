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


## ======================
## ---- subsetByTreatment
## ----------------------
