## Navigating this file:
## - Slot section names start with ----
## - Method section names start with ==
## 
## As a result, you can use Ctrl + f to find the slot or method you are looking
## for quickly, assuming you know its name.
## 
## For example Ctrl + f '== molecularProfiles' would take you the molecularProfiles
## method, while Ctrl +f '---- molecularProfiles' would take you to the slot
## section.

#' @include PharmacoSet-class.R allGenerics.R
NULL

## Variables for dynamic inheritted roxygen2 docs

.local_class <- 'PharmacoSet'
.local_data <- 'CCLEsmall'

#' @title .parseToRoxygen
#'
#' Helper for metaprogramming roxygen2 documentation
#'
#' @description Takes a string block of roxygen2 tags sepearated by new-line 
#'   characteres and parses it to the appropriate format for the @eval tag,
#'   subtituting any string in { } for the argument of the same name in `...`.
#'
#' @keywords internal
#' @export
#' @noRd
.parseToRoxygen <- function(string, ...) {
    unlist(strsplit(
        with(list(...), glue::glue(string)),
    '\n'))
}

# =======================================
# Accessor Method Documentation Object
# ---------------------------------------


#' @name PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_accessors(class_=.local_class)
#' @eval .parseToRoxygen(
#'      "@examples data({data_})
#'      ", data_=.local_data)
#' @importFrom methods callNextMethod
NULL


# ======================================
# Accessor Methods
# --------------------------------------



## ==============
## ---- drug slot
## --------------


##
## == drugInfo

#' @noRd
.docs_PharmacoSet_get_drugInfo <- function(...) .parseToRoxygen(
    "
    @details
    
    ## drug slot accessors
    
    __drugInfo__: `data.frame` Retrieve the drug metadata from a 
    `{class_}` objects `@drug` slot.

    @examples

    ## drug slot

    drugInfo({data_})
    
    @md
    @aliases drugInfo{class_}-method
    @exportMethod drugInfo
    ",
    ...
)

#' @rdname PharmacoSet-accessors
#' @eval .docs_PharmacoSet_get_drugInfo(class_=.local_class, data_=.local_data)
setMethod(drugInfo, signature=c(object='PharmacoSet'), function(object) {
    object@drug
})

#' @noRd
.docs_PharmacoSet_set_drugInfo <- function(...) .parseToRoxygen(
    "
    @details

    __drugInfo__: Update the `@drug` slot of a `{class_}` object.
    - value: `data.frame` of updated drug metadata to assign to a 
    `{class_}` objects `@drug` slot.

    @examples
    drugInfo({data_}) <- drugInfo({data_})

    @md
    @aliases drugInfo<-{class_},data.frame-method
    @exportMethod drugInfo<-
    ",
    ...
)

#' @rdname PharmacoSet-accessors
#' @eval .docs_PharmacoSet_set_drugInfo(class_=.local_class, data_=.local_data)
setReplaceMethod("drugInfo", signature(object="PharmacoSet", 
    value="data.frame"), function(object, value) 
{
    object@drug <- value
    object
})


##
## == drugNames


#' @noRd
.docs_PharmacoSet_get_drugNames <- function(...) .parseToRoxygen(
    "
    @details
    __drugNames__: `character()` The names of all drugs available in a specified 
    `{class_}` object.

    @examples
    drugNames({data_})

    @md
    @aliases drugNames,{class_}-method
    @exportMethod drugNames
    ",
    ...
)

#' @rdname PharmacoSet-accessors
#' @eval .docs_PharmacoSet_get_drugNames(class_=.local_class, data_=.local_data)
setMethod(drugNames, signature="PharmacoSet", function(object) {
  rownames(drugInfo(object))
})

.docs_PharmacoSet_set_drugNames <- function(...) .parseToRoxygen(
    "
    @details
    __drugNames<-__:  Set the drug names available in a PharmacoSet object.
    - value: A `character` vector of the new drug names. Must have the same
    length equal to `nrow(drugInfo(object))`.

    @examples
    drugNames({data_}) <- drugNames({data_})

    @md
    @aliases drugNames<-,{class_},character-method
    @exportMethod drugNames<-
    ",
    ...
)

#' @rdname PharmacoSet-accessors
#' @eval .docs_PharmacoSet_set_drugNames(class_=.local_class, data_=.local_data)
setReplaceMethod("drugNames", signature(object="PharmacoSet", 
    value="character"), function(object, value) 
{
    object <- updateDrugId(object, value)
    return(object)
})


#### CoreGx inherited methods
####
#### Note: The raw documentation lives in CoreGx, see the functions called
#### in @eval tags for the content of the metaprogrammed roxygen2 docs.
####
#### See .parseToRoxygen method in utils-messages.R file of CoreGx to 
#### create similar metaprogrammed docs.
####
#### Warning: for dynamic docs to work, you must set
#### Roxygen: list(markdown = TRUE, r6=FALSE)
#### in the DESCRPTION file!



## ====================
## ---- annotation slot
## --------------------


##
## == annotation


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_annotation(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx annotation
setMethod('annotation', signature("PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors 
#' @eval CoreGx:::.docs_CoreSet_set_annotation(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx annotation<-
setReplaceMethod("annotation", signature("PharmacoSet", "list"),
    function(object, value) 
{
    callNextMethod(object=object, value=value)
})


##
## == dateCreated


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_dateCreated(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx dateCreated
setMethod('dateCreated', signature("PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_dateCreated(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx dateCreated<-
setReplaceMethod('dateCreated', signature(object="PharmacoSet", value="character"), 
    function(object, value) 
{
    callNextMethod(object=object, value=value)
})


##
## === name


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_name(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx name
setMethod('name', signature("PharmacoSet"), function(object){
    return(object@annotation$name)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_name(class_=.local_class, data_=.local_data)
#' @importMethodsFrom CoreGx name<-
setReplaceMethod('name', signature("PharmacoSet"), function(object, value){
    object@annotation$name <- value
    return(object)
})

## ==============
## ---- cell slot
## --------------


##
## == cellInfo


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_cellInfo(class_=.local_class)
#' @importFrom CoreGx cellInfo
setMethod(cellInfo, "PharmacoSet", function(object){
    callNextMethod(object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_cellInfo(class_=.local_class, 
#' data_=.local_data)
#' @importFrom CoreGx cellInfo<-
setReplaceMethod("cellInfo", signature(object="PharmacoSet",
    value="data.frame"), function(object, value)
{
    callNextMethod(object, value=value)
})


##
## == cellNames


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_cellNames(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx cellNames
setMethod(cellNames, signature("PharmacoSet"), function(object){
    callNextMethod(object)
})


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_cellNames(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx cellNames<-
setReplaceMethod("cellNames", signature(object="PharmacoSet",value="character"), 
  function(object, value)
{
    callNextMethod(object=object, value=value)
})



## ------------------
## ---- curation slot


##
## == curation


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_curation(class_=.local_class, 
#' data_=.local_data, details_="Contains three `data.frame`s, 'cell' with 
#' cell-line ids and 'tissue' with tissue ids and 'drug' with drug ids.")
#' @importMethodsFrom CoreGx curation
setMethod('curation', signature(object="PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_curation(class_=.local_class, 
#' data_=.local_data, details_="For a `PharmacoSet` object the slot should 
#' contain tissue, cell-line and drug id `data.frame`s.")
#' @importMethodsFrom CoreGx curation<-
setReplaceMethod("curation", signature(object="PharmacoSet", value="list"), 
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


## ----------------------
## ---- datasetType slot


#
# == datasetType


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_datasetType(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx datasetType
setMethod("datasetType", signature("PharmacoSet"), function(object) {
    callNextMethod(object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_datasetType(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx datasetType<-
setReplaceMethod("datasetType", signature(object="PharmacoSet", 
    value='character'), function(object, value) 
{
    callNextMethod(object=object, value=value)
})


## ---------------------------
## ---- molecularProfiles slot


##
## == molecularProfiles


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_molecularProfiles(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfiles
setMethod(molecularProfiles, "PharmacoSet", function(object, mDataType, assay) 
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_molecularProfiles(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfiles<-
setReplaceMethod("molecularProfiles", signature(object="PharmacoSet", 
    mDataType ="character", assay="character", value="matrix"),
    function(object, mDataType, assay, value)
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})
setReplaceMethod("molecularProfiles",
    signature(object="PharmacoSet", mDataType ="character", assay="missing", 
        value="matrix"), function(object, mDataType, assay, value) 
{
    callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})


##
## == featureInfo


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_featureInfo(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx featureInfo
setMethod(featureInfo, "PharmacoSet", function(object, mDataType) {
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_featureInfo(class_=.local_class, 
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx featureInfo<-
setReplaceMethod("featureInfo", signature(object="PharmacoSet", 
    mDataType ="character",value="data.frame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})
setReplaceMethod("featureInfo", signature(object="PharmacoSet", 
    mDataType ="character",value="DataFrame"), 
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})



##
## == phenoInfo


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_phenoInfo(class_=.local_class, 
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx phenoInfo
setMethod('phenoInfo', signature(object='PharmacoSet', mDataType='character'), 
    function(object, mDataType)
{
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_phenoInfo(class_=.local_class, 
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx phenoInfo<-
setReplaceMethod("phenoInfo", signature(object="PharmacoSet", 
    mDataType ="character", value="data.frame"), 
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})
setReplaceMethod("phenoInfo", signature(object="PharmacoSet", 
    mDataType ="character", value="DataFrame"),
    function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})


##
## == fNames


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_fNames(class_=.local_class, 
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx fNames
setMethod('fNames', signature(object='PharmacoSet', mDataType='character'), 
    function(object, mDataType)
{
    callNextMethod(object=object, mDataType=mDataType)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_fNames(class_=.local_class, 
#' data_=.local_data, mDataType_='rna')
#' @importMethodsFrom CoreGx fNames<-
setReplaceMethod('fNames', signature(object='PharmacoSet', mDataType='character',
    value='character'), function(object, mDataType, value)
{
    callNextMethod(object=object, mDataType=mDataType, value=value)
})


##
## == mDataNames


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_mDataNames(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx mDataNames
setMethod("mDataNames", "PharmacoSet", function(object){
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_mDataNames(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx mDataNames<-
setReplaceMethod("mDataNames", "PharmacoSet", function(object, value){
    callNextMethod(object=object, value=value)
})



##
## == molecularProfilesSlot


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_molecularProfilesSlot(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfilesSlot
setMethod("molecularProfilesSlot", signature("PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setClassUnion('list_or_MAE', c('list', 'MultiAssayExperiment'))

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_molecularProfilesSlot(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx molecularProfilesSlot<-
setReplaceMethod("molecularProfilesSlot", signature("PharmacoSet", "list_or_MAE"),
    function(object, value) 
{
    callNextMethod(object=object, value=value)
})


# ---------------------
## ---- sensitivity slot


##
## == sensitivityInfo

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityInfo(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityInfo
setMethod('sensitivityInfo', signature("PharmacoSet"),
    function(object, dimension, ...) 
{
    callNextMethod(object=object, dimension=dimension, ...)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityInfo(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityInfo<-
setReplaceMethod("sensitivityInfo", signature(object="PharmacoSet", 
    value="data.frame"), function(object, dimension, ..., value) 
{
    callNextMethod(object=object, dimension=dimension, ..., value=value)
})


##
## == sensitvityMeasures


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityMeasures(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityMeasures
setMethod('sensitivityMeasures', signature(object="PharmacoSet"), 
    function(object)
{
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitityMeasures(class_=.local_class, 
#' data_=.local_data) 
setReplaceMethod('sensitivityMeasures',
    signature(object='PharmacoSet', value='character'), function(object, value) 
{
    callNextMethod(object=object, value=value)
})


##
## == sensitivityProfiles


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityProfiles(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityProfiles
setMethod('sensitivityProfiles', signature(object="PharmacoSet"), function(object) 
{
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityProfiles(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityProfiles<-
setReplaceMethod("sensitivityProfiles",
    signature(object="PharmacoSet", value="data.frame"),
    function(object, value) 
{
    callNextMethod(object=object, value=value)
})


#
# == sensitivityRaw


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivityRaw(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityRaw
setMethod("sensitivityRaw", signature("PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensitivityRaw(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensitivityRaw<-
setReplaceMethod('sensitivityRaw', signature("PharmacoSet", "array"), 
    function(object, value) 
{
    callNextMethod(object=object, value=value)
})


#
# == sensitivitySlot


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensitivitySlot(class_=.local_class, 
#'   data_=.local_data)
#' @importMethodsFrom CoreGx sensitivitySlot
setMethod("sensitivitySlot", signature("PharmacoSet"), function(object) {
    callNextMethod(object=object)
})

#' @importClassesFrom CoreGx LongTable
setClassUnion('list_or_LongTable', c('list', 'LongTable'))

#' @rdname PharmacoSet-accessors
#' @importMethodsFrom CoreGx sensitivitySlot<-
#' @eval CoreGx:::.docs_CoreSet_set_sensitivitySlot(class_=.local_class, 
#' data_=.local_data)
setReplaceMethod('sensitivitySlot', signature(object='PharmacoSet', 
    value='list_or_LongTable'), function(object, value) 
{
    callNextMethod(object=object, value=value)
})


##
## == sensNumber


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_sensNumber(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensNumber
setMethod('sensNumber', "PharmacoSet", function(object){
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_sensNumber(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx sensNumber<-
setReplaceMethod('sensNumber', signature(object="PharmacoSet", value="matrix"), 
    function(object, value)
{
    callNextMethod(object=object, value=value)
})


## ======================
## ---- perturbation slot


##
## == pertNumber 


#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_get_pertNumber(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx pertNumber
setMethod('pertNumber', signature(object='PharmacoSet'), function(object) {
    callNextMethod(object=object)
})

#' @rdname PharmacoSet-accessors
#' @eval CoreGx:::.docs_CoreSet_set_pertNumber(class_=.local_class, 
#' data_=.local_data)
#' @importMethodsFrom CoreGx pertNumber<-
setReplaceMethod('pertNumber', signature(object='PharmacoSet', value="array"), 
    function(object, value)
{
    callNextMethod(object=object, value=value)
})