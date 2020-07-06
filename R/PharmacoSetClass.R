#' A Class to Contain PharmacoGenomic datasets together with their curations
#' 
#' The PharmacoSet (pSet) class was developed to contain and organise large 
#' PharmacoGenomic datasets, and aid in their metanalysis. It was designed 
#' primarily to allow bioinformaticians and biologists to work with data at the 
#' level of genes, drugs and cell lines, providing a more naturally intuitive 
#' interface and simplifying analyses between several datasets. As such, it was 
#' designed to be flexible enough to hold datasets of two different natures 
#' while providing a common interface. The class can accomidate datasets 
#' containing both drug dose response data, as well as datasets contaning 
#' genetic profiles of cell lines pre and post treatement with compounds, known 
#' respecitively as sensitivity and perturbation datasets.
#' 
#' @param object A \code{PharmacoSet} object
#' @param mDataType A \code{character} with the type of molecular data to 
#'   return/update
#' @param value A replacement value
#' 
#' @slot annotation A \code{list} of annotation data about the PharmacoSet,
#'    including the \code{$name} and the session information for how the object
#'    was creating, detailing the exact versions of R and all the packages used
#' @slot molecularProfiles A \code{list} containing \code{SummarizedExperiment} 
#'   type object for holding data for RNA, DNA, SNP and CNV 
#'   measurements, with associated \code{fData} and \code{pData} 
#'   containing the row and column metadata
#' @slot cell A \code{data.frame} containing the annotations for all the cell 
#'   lines profiled in the data set, across all data types
#' @slot drug A \code{data.frame} containg the annotations for all the drugs 
#'   profiled in the data set, across all data types
#' @slot sensitivity A \code{list} containing all the data for the sensitivity 
#'   experiments, including \code{$info}, a \code{data.frame} containing the 
#'   experimental info,\code{$raw} a 3D \code{array} containing raw data, 
#'   \code{$profiles}, a \code{data.frame} containing sensitivity profiles 
#'   statistics, and \code{$n}, a \code{data.frame} detailing the number of 
#'   experiments for each cell-drug pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame} 
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for \code{$drug}, 
#'   \code{cell}, \code{tissue} names  used in the data set to universal 
#'   identifiers used between different PharmacoSet objects
#' @slot datasetType A \code{character} string of 'sensitivity', 
#'   'perturbation', or both detailing what type of data can be found in the 
#'   PharmacoSet, for proper processing of the data
#'  
#' @importClassesFrom CoreGx CoreSet
#' @return An object of the PharmacoSet class
.PharmacoSet <- setClass("PharmacoSet", slots = list(drug="data.frame"), 
                         contains="CoreSet")


# The default constructor above does a poor job of explaining the required 
# structure of a PharmacoSet. The constructor function defined below guides the 
# user into providing the required components of the curation and senstivity 
# lists and hides the annotation slot which the user does not need to manually 
# fill. This also follows the design of the Expression Set class.

#####
# CONSTRUCTOR -----
#####

#' PharmacoSet constructor
#' 
#' A constructor that simplifies the process of creating PharmacoSets, as well 
#' as creates empty objects for data not provided to the constructor. Only
#' objects returned by this constructor are expected to work with the PharmacoSet
#' methods. For a much more detailed instruction on creating PharmacoSets, please
#' see the "CreatingPharmacoSet" vignette.
#' 
#' @examples  
#' ## For help creating a PharmacoSet object, please see the following vignette:
#' browseVignettes("PharmacoGx")
#' 
##TODO:: Determine how to generalise the constructor documentation in CoreGx
## to make sense with all three packages which depend on it. For now it says
## CoreSet instead of PharmacoSet
##TODO:: Determine if there is any way to execture R code when making roxygen2
## documentation. Then we could use a variable to fill in the class for each
## package.
#' @inheritParams CoreGx::CoreSet
# @param name A \code{character} string detailing the name of the dataset
# @param molecularProfiles A \code{list} of SummarizedExperiment objects containing
#   molecular profiles for each data type.
# @param cell A \code{data.frame} containg the annotations for all the cell
#   lines profiled in the data set, across all data types
#' @param drug A \code{data.frame} containg the annotations for all the drugs
#   profiled in the data set, across all data types
# @param sensitivityInfo A \code{data.frame} containing the information for the
#   sensitivity experiments
# @param sensitivityRaw A 3 Dimensional \code{array} contaning the raw drug
#   dose response data for the sensitivity experiments
# @param sensitivityProfiles \code{data.frame} containing drug sensitivity profile
#   statistics such as IC50 and AUC
# @param sensitivityN,perturbationN A \code{data.frame} summarizing the
#   available sensitivity/perturbation data
#' @param curationDrug,curationCell,curationTissue A \code{data.frame} mapping
#'   the names for drugs, cells and tissues used in the data set to universal
#'   identifiers used between different PharmacoSet objects
# @param datasetType A \code{character} string of 'sensitivity',
#   'preturbation', or both detailing what type of data can be found in the
#   PharmacoSet, for proper processing of the data
# @param verify \code{boolean} Should the function verify the PharmacoSet and
#   print out any errors it finds after construction?
#' 
#' @return An object of class PharmacoSet
#
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment rowData colData assay assays assayNames Assays
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom CoreGx CoreSet
#' 
#' @export
PharmacoSet <-  function(name,
                          molecularProfiles=list(), 
                          cell=data.frame(), 
                          drug=data.frame(),
                          sensitivityInfo=data.frame(),
                          sensitivityRaw=array(dim=c(0,0,0)), 
                          sensitivityProfiles=matrix(), 
                          sensitivityN=matrix(nrow=0, ncol=0), 
                          perturbationN=array(NA, dim=c(0,0,0)), 
                          curationDrug=data.frame(), 
                          curationCell = data.frame(), 
                          curationTissue = data.frame(), 
                          datasetType=c("sensitivity", "perturbation", "both"),
                          verify = TRUE
                         )
{
    datasetType <- match.arg(datasetType)
    
    annotation <- list()
    annotation$name <- as.character(name)
    annotation$dateCreated <- date()
    annotation$sessionInfo <- sessionInfo()
    annotation$call <- match.call()
    
    ## TODO:: If the colnames and rownames are not found below, it will fill with NAs. This is undersirable behaviour.
    #molecularProfiles <- list("dna"=dna, "rna"=rna, "snp"=snp, "cnv"=cnv)
    ## TODO:: Determine if I should use SummarizedExperiment construtor here?
    for (i in seq_along(molecularProfiles)){
        if (!is(molecularProfiles[[i]], "SummarizedExperiment")) {
            stop(sprintf("Please provide the %s data as a SummarizedExperiment", 
                         names(molecularProfiles[i])))
        }else{
      rowData(molecularProfiles[[i]]) <- 
        rowData(molecularProfiles[[i]])[rownames(assays(molecularProfiles[[i]])[[1]]), , drop=FALSE]
      colData(molecularProfiles[[i]]) <- 
        colData(molecularProfiles[[i]])[colnames(assays(molecularProfiles[[i]])[[1]]), , drop=FALSE]
        }
    
    }
    #if (!is(cell, "data.frame")) {
    #    stop("Please provide the cell line annotations as a data frame.")
    #}
    #if (!is(drug, "data.frame")) {
    #    stop("Please provide the drug annotations as a data frame.")
    #}
    
    sensitivity <- list()
    
    if (!all(rownames(sensitivityInfo) == rownames(sensitivityProfiles) & 
             rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])){
        stop("Please ensure all the row names match between the sensitivity data.")
    }
    
    sensitivity$info <- as.data.frame(sensitivityInfo, stringsAsFactors = FALSE)
    sensitivity$raw <- sensitivityRaw
    sensitivity$profiles <- as.data.frame(sensitivityProfiles, stringsAsFactors = FALSE)
    sensitivity$n <- sensitivityN
    
    curation <- list()
    curation$drug <- as.data.frame(curationDrug, stringsAsFactors = FALSE)
    curation$cell <- as.data.frame(curationCell, stringsAsFactors = FALSE)
    curation$tissue <- as.data.frame(curationTissue, stringsAsFactors = FALSE)
    ### TODO:: Make sure to fix the curation to check for matching row names to the drug and cell line matrices!!!!!!
    
    perturbation <- list()
    perturbation$n <- perturbationN
    if (datasetType == "perturbation" || datasetType == "both") {
        perturbation$info <- "The metadata for the perturbation experiments is 
          available for each molecular type by calling the appropriate info 
          function. \n For example, for RNA transcriptome perturbations, 
          the metadata can be accessed using rnaInfo(object)."
    } else {
        perturbation$info <- "Not a perturbation dataset."
    }
    
    pSet  <- .PharmacoSet(annotation=annotation, molecularProfiles=molecularProfiles, cell=as.data.frame(cell), drug=as.data.frame(drug), datasetType=datasetType, sensitivity=sensitivity, perturbation=perturbation, curation=curation)
    if (verify) { checkPsetStructure(pSet) }
  if(length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
    pSet@sensitivity$n <- .summarizeSensitivityNumbers(pSet)
  }
    if(length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
      pSet@perturbation$n <- .summarizePerturbationNumbers(pSet)
    }
  return(pSet)
}

##TODO:: Figure out how to properly inherit params from CoreGx

#####
# CELL SLOT GETTERS/SETTERS ----
#####

#' cellInfo
#'
#' Return cell line metadata from a object
#'
#' @examples
#' data(CCLEsmall)
#' cellInf <- cellInfo(CCLEsmall)
#'
#' @param object The \code{PharmacoSet} to retrieve cell info from
#'
#' @return a \code{data.frame} with the cell annotations
#' 
#' @describeIn PharmacoSet
#' 
#' @importFrom CoreGx cellInfo
#' @importFrom methods callNextMethod
#' @export
setMethod(cellInfo, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

#' cellInfo<- Generic
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
#' @export
setReplaceMethod("cellInfo", signature = signature(object="PharmacoSet", value="data.frame"), function(object, value){
  callNextMethod(object=object, value=value)
})

#####
# DRUG SLOT GETTERS/SETTERS ----
#####

#' drugInfo Generic
#' 
#' Generic for drugInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' drugInf <- drugInfo(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to retrieve drug info from
#' 
#' @return a \code{data.frame} with the drug annotations
setGeneric("drugInfo", function(object) standardGeneric("drugInfo"))
#'
#' @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' 
#' @export
setMethod(drugInfo, "PharmacoSet", function(object){
  object@drug
})

#' drugInfo<- Generic
#' 
#' Generic for drugInfo replace method
#' 
#' @examples
#' data(CCLEsmall)
#' drugInfo(CCLEsmall) <- drugInfo(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to replace drug info in
#' @param value A \code{data.frame} with the new drug annotations
#' 
#' @return Updated \code{PharmacoSet}
#' 
setGeneric('drugInfo<-', function(object, value) standardGeneric('drugInfo<-'))
#' @describeIn PharmacoSet Update the drug annotations
#' 
#' @export
setReplaceMethod("drugInfo", signature = signature(object="PharmacoSet", value="data.frame"), function(object, value) {
  object@drug <- value
  object
})

#####
# MOLECULAR PROFILES SLOT GETTERS/SETTERS ----
#####

#' phenoInfo Generic
#' 
#' Generic for phenoInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' phenoInf <- phenoInfo(CCLEsmall, mDataType="rna")
#' 
#' @param object The \code{PharmacoSet} to retrieve rna annotations from
#' @param mDataType the type of molecular data 
#' 
#' @return a \code{data.frame} with the experiment info
#' 
#' @describeIn PharmacoSet Return the experiment info from the given type of 
#'   molecular data in PharmacoSet
#'   
#' @importFrom CoreGx phenoInfo
#' @importFrom methods callNextMethod
#'
#' @export
setMethod(phenoInfo, "PharmacoSet", function(object, mDataType){
    callNextMethod(object=object, mDataType=mDataType)
})

#' phenoInfo<- Generic
#' 
#' Generic for phenoInfo replace method
#' 
#' @examples
#' data(CCLEsmall)
#' phenoInfo(CCLEsmall, mDataType="rna") <- phenoInfo(CCLEsmall, mDataType="rna")
#' 
#' @param object The \code{PharmacoSet} to retrieve molecular experiment annotations from
#' @param mDataType the type of molecular data 
#' @param value a \code{dataframe}  with the new experiment annotations
#' 
#' @return The updated \code{PharmacoSet}
#'
#' @describeIn PharmacoSet Update the given type of molecular data experiment info in the PharmacoSet
#' 
#' @importFrom methods callNextMethod
#' @importFrom CoreGx phenoInfo<-
#' @export 
setReplaceMethod("phenoInfo", signature = signature(object="PharmacoSet", mDataType ="character", value="DataFrame"), function(object, mDataType, value){
  callNextMethod(object=object, mDataType=mDataType, value=value)
})

#' molecularProfiles Generic
#' 
#' Generic for molecularProfiles method 
#' 
#' @examples
#' data(CCLEsmall)
#' molProf <- molecularProfiles(CCLEsmall, "rna")
#' 
#' @param object The \code{PharmacoSet} to retrieve molecular profiles from
#' @param mDataType \code{character} The type of molecular data
#' @param assay \code{character} Name of the desired assay; if excluded defaults to first assay
#'   in the SummarizedExperiment for the given mDataType
#' 
#' @return a \code{matrix} of data for the given mDataType and assay
#' 
#' @describeIn PharmacoSet Return the given type of molecular data from the PharmacoSet 
#' 
#' @importFrom CoreGx molecularProfiles
#' @importFrom methods callNextMethod
#' @export
setMethod(molecularProfiles, "PharmacoSet", function(object, mDataType, assay){
  callNextMethod(object=object, mDataType=mDataType, assay=assay)
})

#' molecularProfiles<- Generic
#' 
#' Generic for molecularProfiles replace method
#' 
#' @examples
#' data(CCLEsmall)
#' molecularProfiles(CCLEsmall, "rna") <- molecularProfiles(CCLEsmall, "rna")
#' 
#' @param object The \code{PharmacoSet} to replace molecular profiles in
#' @param mDataType The type of molecular data to be updated
#' @param assay \code{character} Name or index of the assay data to return
#' @param value A \code{matrix} with the new profiles
#' 
#' @return Updated \code{PharmacoSet}
#' 
#' @describeIn PharmacoSet Update the given type of molecular data from the PharmacoSet 
#' 
#' @importFrom methods callNextMethod
#' @importFrom CoreGx molecularProfiles<-
#' @export
setReplaceMethod("molecularProfiles", signature = signature(object="PharmacoSet", mDataType ="character", assay="character", value="matrix"), function(object, mDataType, assay, value){
  callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})
#' @describeIn PharmacoSet Update the given type of molecular data from the PharmacoSet 
#' @export
setReplaceMethod("molecularProfiles", signature = signature(object="PharmacoSet", mDataType ="character", assay="missing", value="matrix"), function(object, mDataType, assay, value){
  callNextMethod(object=object, mDataType=mDataType, assay=assay, value=value)
})

#' featureInfo Generic
#' 
#' Generic for featureInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' featInf <- featureInfo(CCLEsmall, "rna")
#' 
#' @param object The \code{PharmacoSet} to retrieve feature annotations from
#' @param mDataType the type of molecular data 
#' @return a \code{data.frame} with the experiment info
#'
#' @describeIn PharmacoSet Return the feature info for the given molecular data
#' 
#' @importFrom methods callNextMethod
#' @importFrom CoreGx featureInfo
#' @export
setMethod(featureInfo, "PharmacoSet", function(object, mDataType){
  callNextMethod(object=object, mDataType=mDataType)
})

#' featureInfo<- Generic
#' 
#' Generic for featureInfo replace method
#' 
#' @examples
#' data(CCLEsmall)
#' featureInfo(CCLEsmall, "rna") <- featureInfo(CCLEsmall, "rna")
#' 
#' @param object The \code{PharmacoSet} to replace gene annotations in
#' @param mDataType The type of molecular data to be updated
#' @param value A \code{DataFrame} with the new feature annotations
#' 
#' @return Updated \code{PharmacoSet}
#' 
#' @importFrom CoreGx featureInfo<-
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Replace the gene info for the molecular data
#' 
#' @export
setReplaceMethod("featureInfo", signature = signature(object="PharmacoSet", mDataType ="character", value="DataFrame"), function(object, mDataType, value){
  callNextMethod(object=object, mDataType=mDataType, value=value)
})

#####
# SENSITIVITY SLOT GETTERS/SETTERS
#####

#' sensitivityInfo Generic
#' 
#' Get the senstivity information DataFrame from a object object
#' 
#' @examples
#' data(CCLEsmall)
#' sensInf <- sensitivityInfo(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to retrieve sensitivity experiment annotations from
#' 
#' @return a \code{DataFrame} with the experiment info
#'
#' @importFrom CoreGx sensitivityInfo
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the drug dose sensitivity experiment info
#'
#' @export
setMethod(sensitivityInfo, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

#' sensitivityInfo<- Generic
#' 
#' Set the sensitivityInfo DataFrame in a object object
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{DataFrame} with the new sensitivity annotations
#' 
#' @return Updated \code{PharmacoSet} 
#' 
#' @importFrom CoreGx sensitivityInfo<-
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Update the sensitivity experiment info
#' @export
setReplaceMethod("sensitivityInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){
  callNextMethod(object=object, value=value)
})


#' sensitivityProfiles
#' 
#' Get the sensitivityProfiles data.frame from a object object
#' 
#' @examples
#' data(CCLEsmall)
#' sensProf <- sensitivityProfiles(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to retrieve sensitivity experiment data from
#' 
#' @return a \code{data.frame} with the experiment info
#'
#' @importFrom CoreGx sensitivityProfiles
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Return the phenotypic data for the drug dose sensitivity
#' @export
setMethod(sensitivityProfiles, "PharmacoSet", function(object) {
  callNextMethod(object=object)
})

#' sensitivityProfiles<-
#' 
#' Set the sensitivityProfiles in a object object
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityProfiles(CCLEsmall) <- sensitivityProfiles(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{data.frame} with the new sensitivity profiles. If a matrix object is passed in, converted to data.frame before assignment
#' 
#' @return Updated \code{PharmacoSet} 
#' 
#' @importFrom CoreGx sensitivityProfiles<-
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet", value="data.frame"), function(object, value){
  callNextMethod(object=object, value=value)
})
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet", value="matrix"), function(object, value){
  callNextMethod(object=object, value=value)
})

#' sensitivityMeasures Generic
#' 
#' Get the types of sensitivity measurements from a object object
#' 
#' @examples
#' data(CCLEsmall)
#' sensMeas <- sensitivityMeasures(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} 
#' 
#' @return A \code{character} vector of all the available sensitivity measures
#'
#' @importFrom CoreGx sensitivityMeasures
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#' @export
setMethod(sensitivityMeasures, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

##TODO:: Arrange slot accessors to be together!

#' drugNames Generic
#' 
#' A generic for the drugNames method
#' 
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to return drug names from
#' @return A vector of the drug names used in the PharmacoSet
setGeneric("drugNames", function(object) standardGeneric("drugNames"))
#' @describeIn PharmacoSet Return the names of the drugs used in the PharmacoSet
#' @export
setMethod(drugNames, "PharmacoSet", function(object){
  rownames(drugInfo(object))
})

#' drugNames<- Generic
#' 
#' A generic for the drugNames replacement method
#' 
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall) <- drugNames(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{character} vector of the new drug names
#' @return Updated \code{PharmacoSet} 
setGeneric("drugNames<-", function(object, value) standardGeneric("drugNames<-"))
#' @describeIn PharmacoSet Update the drug names used in the dataset
#' @export
setReplaceMethod("drugNames", signature = signature(object="PharmacoSet",value="character"), function(object, value){
    object <- updateDrugId(object, value)
    return(object)
    })

#' cellNames Generic
#' 
#' A generic for the cellNames method
#' 
#' @examples
#' data(CCLEsmall)
#' cellNames(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to return cell names from
#' @return A vector of the cell names used in the PharmacoSet
#' 
#' @importFrom CoreGx cellNames
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Return the cell names used in the dataset
#' @export
#' 
setMethod(cellNames, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

#' cellNames<- Generic
#' 
#' A generic for the cellNames replacement method
#' 
#' @examples
#' data(CCLEsmall)
#' cellNames(CCLEsmall) <- cellNames(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{character} vector of the new cell names
#' 
#' @return Updated \code{PharmacoSet}
#'  
#' @importFrom CoreGx cellNames<-
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Update the cell names used in the dataset
#' @export
setReplaceMethod("cellNames", signature = signature(object="PharmacoSet",value="character"), function(object, value){
    callNextMethod(object=object, value=value)
})

#' fNames
#' 
#' Return the feature names for the specified molecular data type
#' 
#' @examples
#' data(CCLEsmall)
#' fNames(CCLEsmall, "rna")
#'
#' @param object The \code{PharmacoSet} 
#' @param mDataType The molecular data type to return feature names for
#' @return A \code{character} vector of the feature names
#' 
#' @importFrom CoreGx fNames
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Return the feature names used in the dataset
#' @export
setMethod(fNames, signature=signature(object="PharmacoSet", mDataType="character"), function(object, mDataType){
  callNextMethod(object=object, mDataType=mDataType)
})

#' fNames<-
#'
#' Setter for the feature names of a \code{SummarizedExperiment} in the 
#'   molecularProfiles slot
#'
#' @examples
#' data(CCLEsmall)
#' fNames(CCLEsmall, 'rna') <- fNames(CCLEsmall, 'rna')
#' 
#' @param object The \code{PharmacoSet} object to update
#' @param mDataType The molecular data type to update
#' @param value A \code{character} vector of the new cell names
#' @return Updated \code{PharmacoSet} 
#'
#' @importFrom CoreGx fNames<-
#' @importFrom methods callNextMethod
#'
#' @export
setReplaceMethod("fNames", 
                 signature = signature(object="PharmacoSet",
                                       mDataType='character', 
                                       value="character"), 
                 function(object, mDataType, value)
{
  callNextMethod(object=object, mDataType=mDataType, value=value)
})

#' dateCreated Generic
#' 
#' A generic for the dateCreated method
#' 
#' @examples
#' data(CCLEsmall)
#' dateCreated(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} 
#' @return The date the PharmacoSet was created
#'
#' @importFrom CoreGx dateCreated
#' @importFrom methods callNextMethod
#'
#' @describeIn PharmacoSet Return the date the PharmacoSet was created
#' @export
setMethod(dateCreated, "PharmacoSet", function(object) {
  callNextMethod(object=object)
})

#' name getter method
#' 
#' Retrun the name of the PharmacoSet object
#' 
#' @examples
#' data(CCLEsmall)
#' name(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet}
#' @return The name of the PharmacoSet
#'
#' @importFrom CoreGx name
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Return the name of the PharmacoSet 
#' @export
setMethod(name, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

#' pertNumber Generic
#' 
#' A generic for the pertNumber method
#' 
#' @examples
#' data(CCLEsmall)
#' pertNumber(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} 
#' @return A 3D \code{array} with the number of perturbation experiments per 
#'   drug and cell line, and data type
#'
#' @importFrom CoreGx pertNumber
#' @importFrom methods callNextMethod
#' 
#' @describeIn PharmacoSet Return the summary of available perturbation
#'   experiments
#' @export
setMethod(pertNumber, "PharmacoSet", function(object){
  callNextMethod(object=object)
})


#' sensNumber Generic
#' 
#' A generic for the sensNumber method
#' 
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} 
#' @return A \code{data.frame} with the number of sensitivity experiments per 
#'   drug and cell line
#'
#' @describeIn PharmacoSet Return the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber
#' @importFrom methods callNextMethod
#' @export
setMethod(sensNumber, "PharmacoSet", function(object){
  callNextMethod(object=object)
})

#' pertNumber<- Generic
#' 
#' A generic for the pertNumber method
#' 
#' @examples
#' data(CCLEsmall)
#' pertNumber(CCLEsmall) <- pertNumber(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} 
#' @param value A new 3D \code{array} with the number of perturbation experiments per drug and cell line, and data type
#' 
#' @return The updated \code{PharmacoSet} 
#'
#' @describeIn PharmacoSet Update the summary of available perturbation
#'   experiments
#'   
#' @importFrom CoreGx pertNumber<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod('pertNumber', signature = signature(object="PharmacoSet", value="array"), function(object, value){
  callNextMethod(object=object, value=value)
})

#' sensNumber<- Generic
#' 
#' A generic for the sensNumber method
#' 
#' 
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall) <- sensNumber(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} 
#' @param value A new \code{data.frame} with the number of sensitivity experiments per drug and cell line
#' 
#' @return The updated \code{PharmacoSet} 
#'
#' @describeIn PharmacoSet Update the summary of available sensitivity
#'   experiments
#'
#' @importFrom CoreGx sensNumber<-
#' @importFrom methods callNextMethod
#' @export
setReplaceMethod('sensNumber', signature = signature(object="PharmacoSet",value="matrix"), function(object, value){
  callNextMethod(object=object, value=value)
})

#' Show a PharamcoSet
#' 
#' @param object \code{PharmacoSet}
#' 
#' @examples
#' data(CCLEsmall)
#' CCLEsmall
#' 
#' @return Prints the PharmacoSet object to the output stream, and returns 
#'   invisible NULL. 
#' 
#'  @importFrom CoreGx show
#'  @importFrom methods callNextMethod
#'
#' @export
setMethod("show", signature=signature(object="PharmacoSet"), function(object) {
  callNextMethod(object=object)
})

#' mDataNames
#' 
#' Returns the molecular data names for the PharmacoSet.
#' 
#' @examples
#' data(CCLEsmall)
#' mDataNames(CCLEsmall)
#' 
#' @param object PharamcoSet object
#' 
#' @return Vector of names of the molecular data types
#' 
#' @importFrom CoreGx mDataNames
#' @importFrom methods callNextMethod
#' @export
setMethod('mDataNames', "PharmacoSet", function(object) {
  callNextMethod(object=object)
})


#'`[`
#'
#'@param x object
#'@param i Cell lines to keep in object
#'@param j Drugs to keep in object
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted object
#'@export
setMethod(`[`, "PharmacoSet", function(x, i, j, ..., drop = FALSE){
  if(is.character(i)&&is.character(j)){
    return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
  } 
  else if(is.numeric(i) && is.numeric(j) && all(as.integer(i)==i) && all(as.integer(j)==j)){
    return(subsetTo(x, cells=cellNames(x)[i], drugs=drugNames(x)[j],  molecular.data.cells=cellNames(x)[i]))
  }
})

#' Get the dimensions of a PharmacoSet
#' 
#' @param x PharmacoSet
#' @return A named vector with the number of Cells and Drugs in the PharmacoSet
#' @export
setMethod("dim", signature=signature(x="PharmacoSet"), function(x){
  return(c(Cells=length(cellNames(x)), Drugs=length(drugNames(x))))
})

## FIXED? TODO:: Subset function breaks if it doesnt find cell line in sensitivity info
#' A function to subset a PharmacoSet to data containing only specified drugs, cells and genes
#' 
#' This is the prefered method of subsetting a PharmacoSet. This function allows
#' abstraction of the data to the level of biologically relevant objects: drugs
#' and cells. The function will automatically go through all of the
#' combined data in the PharmacoSet and ensure only the requested drugs
#' and cell lines are found in any of the slots. This allows quickly picking out
#' all the experiments for a drug or cell of interest, as well removes the need
#' to keep track of all the metadata conventions between different datasets.
#' 
#' @examples
#' data(CCLEsmall)
#' CCLEdrugs  <- drugNames(CCLEsmall)
#' CCLEcells <- cellNames(CCLEsmall)
#' pSet <- subsetTo(CCLEsmall, drugs = CCLEdrugs[1], cells = CCLEcells[1])
#' pSet
#' 
#' @param object A \code{PharmacoSet} to be subsetted
#' @param cells A list or vector of cell names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all cells will be left in
#'   the dataset.
#' @param drugs A list or vector of drug names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all drugs will be left in
#'   the dataset.
#' @param molecular.data.cells A list or vector of cell names to keep in the
#'   molecular data
#' @param keep.controls If the dataset has perturbation type experiments, should
#'   the controls be kept in the dataset? Defaults to true.
#' @param ... Other arguments passed by other function within the package
#' 
#' @return A PharmacoSet with only the selected drugs and cells
#' 
#' @importMethodsFrom CoreGx subsetTo
#' @export
setMethod('subsetTo', signature(object='PharmacoSet'), function(object, cells=NULL, drugs=NULL, molecular.data.cells=NULL,
                     keep.controls=TRUE, ...){
  .subsetToPharmacoSet(object=object, cells=cells, drugs=drugs, molecular.data.cells=molecular.data.cells,
                       keep.controls=keep.controls)
})

#' @importFrom CoreGx .intersectList
#' @keywords internal
.subsetToPharmacoSet <- function(object, cells=NULL, drugs=NULL, molecular.data.cells=NULL,
                     keep.controls=TRUE, ...) {
  drop=FALSE #TODO:: Is this supposed to be here?
  
  adArgs = list(...)
  if ("exps" %in% names(adArgs)) {
  	exps <- adArgs[["exps"]]
  	if(is(exps, "data.frame")) {
  		exps2 <- exps[[name(object)]]
  		names(exps2) <- rownames(exps)
  		exps <- exps2
  	} else{
  		exps <- exps[[name(object)]]
  	}
  }else {
    exps <- NULL
  }
  if(!missing(cells)){
    cells <- unique(cells)
  }
  
  if(!missing(drugs)){
    drugs <- unique(drugs)
  }
  
  if(!missing(molecular.data.cells)){
    molecular.data.cells <- unique(molecular.data.cells)
  }
  
    ### TODO:: implement strict subsetting at this level!!!!
    ### TODO:: refactor this monstrosity of a function into helpers
  
    ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
    ### function in lapply, so it does not recognize them as missing
  
  object@molecularProfiles <- lapply(object@molecularProfiles, function(SE, cells, drugs, molecular.data.cells){
    
    molecular.data.type <- ifelse(length(grep("rna", S4Vectors::metadata(SE)$annotation) > 0), "rna", S4Vectors::metadata(SE)$annotation)
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cells <- molecular.data.cells[[molecular.data.type]]
    }
  
        column_indices <- NULL
  
      if (length(cells)==0 && length(drugs)==0) {
          column_indices <- seq_len(ncol(SE)) # This still returns the number of samples in an SE, but without a label
      }
      if(length(cells)==0 && object@datasetType=="sensitivity") {
        column_indices <- seq_len(ncol(SE))
      }
  
      cell_line_index <- NULL
      if(length(cells)!=0) {
        if (!all(cells %in% cellNames(object))) {
              stop("Some of the cell names passed to function did not match to names in the PharmacoSet. Please ensure you are using cell names as returned by the cellNames function")
        }
          cell_line_index <- which(SummarizedExperiment::colData(SE)[["cellid"]] %in% cells)
        # if (length(na.omit(cell_line_index))==0){
    #       stop("No cell lines matched")
    #     }
      }
      drugs_index <- NULL
      if(object@datasetType=="perturbation" || object@datasetType=="both"){
        if(length(drugs) != 0) {
            if (!all(drugs %in% drugNames(object))){
                  stop("Some of the drug names passed to function did not match to names in the PharmacoSet. Please ensure you are using drug names as returned by the drugNames function")
            }
          drugs_index <- which(SummarizedExperiment::colData(SE)[["drugid"]] %in% drugs)
          # if (length(drugs_index)==0){
    #         stop("No drugs matched")
    #       }
          if(keep.controls) {
            control_indices <- which(SummarizedExperiment::colData(SE)[["xptype"]]=="control")
            drugs_index <- c(drugs_index, control_indices)
          }
        }
      }
  
      if(length(drugs_index) != 0 && length(cell_line_index) != 0) {
        if(length(intersect(drugs_index, cell_line_index)) == 0) {
          stop("This Drug - Cell Line combination was not tested together.")
        }
        column_indices <- intersect(drugs_index, cell_line_index)
      } else {
        if(length(drugs_index) !=0) {
        column_indices <- drugs_index
      }
        if(length(cell_line_index) !=0) {
        column_indices <- cell_line_index
      }
      }
  
      row_indices <- seq_len(nrow(SummarizedExperiment::assay(SE, 1)))
  
      SE <- SE[row_indices, column_indices]
      return(SE)

  }, cells=cells, drugs=drugs, molecular.data.cells=molecular.data.cells)  
  
  if ((object@datasetType == "sensitivity" | object@datasetType == "both") & length(exps) != 0) {
      object@sensitivity$info <- object@sensitivity$info[exps, , drop=drop]
      rownames(object@sensitivity$info) <- names(exps)
      if(length(object@sensitivity$raw) > 0) {
        object@sensitivity$raw <- object@sensitivity$raw[exps, , , drop=drop]
        dimnames(object@sensitivity$raw)[[1]] <- names(exps)
      }
      object@sensitivity$profiles <- object@sensitivity$profiles[exps, , drop=drop]
      rownames(object@sensitivity$profiles) <- names(exps)
      
      object@sensitivity$n <- .summarizeSensitivityNumbers(object)
  }
  else if ((object@datasetType == "sensitivity" | object@datasetType == "both") & (length(drugs) != 0 | length(cells) != 0)) {
    
        drugs_index <- which (sensitivityInfo(object)[, "drugid"] %in% drugs)
        cell_line_index <- which (sensitivityInfo(object)[,"cellid"] %in% cells)
        if (length(drugs_index) !=0 & length(cell_line_index) !=0 ) {
          if (length(intersect(drugs_index, cell_line_index)) == 0) {
            stop("This Drug - Cell Line combination was not tested together.")
          }
          row_indices <- intersect(drugs_index, cell_line_index)
        } else {
          if(length(drugs_index)!=0 & length(cells)==0) {
                row_indices <- drugs_index
          } else {
              if(length(cell_line_index)!=0 & length(drugs)==0){
                  row_indices <- cell_line_index
              } else {
              row_indices <- vector()
              }
          }
       }
        object@sensitivity[names(object@sensitivity)[names(object@sensitivity)!="n"]] <- lapply(object@sensitivity[names(object@sensitivity)[names(object@sensitivity)!="n"]], function(x,i, drop){
            #browser()
          if (length(dim(x))==2){
            return(x[i,,drop=drop])
          }
          if (length(dim(x))==3){
            return(x[i,,,drop=drop])
          }
          }, i=row_indices, drop=drop)
  }
  
	if (length(drugs)==0) {
		if(object@datasetType == "sensitivity" | object@datasetType == "both"){
			drugs <- unique(sensitivityInfo(object)[["drugid"]])
		}
		if(object@datasetType == "perturbation" | object@datasetType == "both"){
			drugs <- union(drugs, na.omit(CoreGx::.unionList(lapply(object@molecularProfiles, function(SE){unique(colData(SE)[["drugid"]])}))))
		}
	}
	if (length(cells)==0) {
		cells <- union(cells, na.omit(CoreGx::.unionList(lapply(object@molecularProfiles, function(SE){unique(colData(SE)[["cellid"]])}))))
        if (object@datasetType =="sensitivity" | object@datasetType == "both"){
            cells <- union(cells, sensitivityInfo(object)[["cellid"]])
        }
	}
	drugInfo(object) <- drugInfo(object)[drugs , , drop=drop]
	cellInfo(object) <- cellInfo(object)[cells , , drop=drop]
	object@curation$drug <- object@curation$drug[drugs , , drop=drop]
	object@curation$cell <- object@curation$cell[cells , , drop=drop]
	object@curation$tissue <- object@curation$tissue[cells , , drop=drop]
	if (object@datasetType == "sensitivity" | object@datasetType == "both"  & length(exps) == 0) {
	  object@sensitivity$n <- object@sensitivity$n[cells, drugs, drop=drop]
	}
	if (object@datasetType == "perturbation" | object@datasetType == "both") {
	  object@perturbation$n <- object@perturbation$n[cells, drugs, , drop=drop]
    }
      return(object)
}
### TODO:: Add updating of sensitivity Number tables
#' @importFrom CoreGx updateCellId
updateCellId <- function(object, new.ids = vector("character")){
  CoreGx::updateCellId(object, new.ids)
}

# updateFeatureNames <- function(object, new.ids = vector("character")){
#
#   if (length(new.ids)!=nrow(cellInfo(object))){
#     stop("Wrong number of cell identifiers")
#   }
#
#   if(object@datasetType=="sensitivity"|object@datasetType=="both"){
#     myx <- match(sensitivityInfo(object)[,"cellid"],rownames(cellInfo(object)))
#     sensitivityInfo(object)[,"cellid"] <- new.ids[myx]
#
#   }
#
#   object@molecularProfiles <- lapply(object@molecularProfiles, function(eset){
#
#     myx <- match(pData(eset)[["cellid"]],rownames(cellInfo(object)))
#     pData(eset)[["cellid"]]  <- new.ids[myx]
#     return(eset)
#       })
#   myx <- match(rownames(object@curation$cell),rownames(cellInfo(object)))
#   rownames(object@curation$cell) <- new.ids[myx]
#   rownames(object@curation$tissue) <- new.ids[myx]
#   if (dim(pertNumber(object))[[1]]>0){
#     myx <- match(dimnames(pertNumber(object))[[1]], rownames(cellInfo(object)))
#     dimnames(pertNumber(object))[[1]] <- new.ids[myx]
#   }
#   if (nrow(sensNumber(object))>0){
#     myx <- match(rownames(sensNumber(object)), rownames(cellInfo(object)))
#     rownames(sensNumber(object)) <- new.ids[myx]
#   }
#   rownames(cellInfo(object)) <- new.ids
#   return(object)
#
# }

### TODO:: Add updating of sensitivity Number tables
updateDrugId <- function(object, new.ids = vector("character")){

  if (length(new.ids)!=nrow(drugInfo(object))){
    stop("Wrong number of drug identifiers")
  }

  if(object@datasetType=="sensitivity"|object@datasetType=="both"){
    myx <- match(sensitivityInfo(object)[,"drugid"],rownames(drugInfo(object)))
    sensitivityInfo(object)[,"drugid"] <- new.ids[myx]

  }
  if(object@datasetType == "perturbation" | object@datasetType == "both"){
    object@molecularProfiles <- lapply(object@molecularProfiles, function(SE) {

      myx <- match(SummarizedExperiment::colData(SE)[["drugid"]],rownames(drugInfo(object)))
      SummarizedExperiment::colData(SE)[["drugid"]]  <- new.ids[myx]
      return(SE)
    })
  }
  

  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateDrugId. Merging old ids into the same identifier")
    
    if(ncol(sensNumber(object))>0){
      sensMatch <- match(colnames(sensNumber(object)), rownames(drugInfo(object)))
    }
    if(dim(pertNumber(object))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(object))[[2]], rownames(drugInfo(object)))
    }
    curMatch <- match(rownames(object@curation$drug),rownames(drugInfo(object)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(object))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(object)[,myx[1]] <- apply(sensNumber(object)[,myx], 1, sum)
        sensNumber(object) <- sensNumber(object)[,-myx[-1]]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(object))[[2]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(object)[,myx[1],] <- apply(pertNumber(object)[,myx,], c(1,3), sum)
        pertNumber(object) <- pertNumber(object)[,-myx[-1],]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      object@curation$drug[myx[1],] <- apply(object@curation$drug[myx,], 2, paste, collapse="///")
      object@curation$drug <- object@curation$drug[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      drugInfo(object)[myx[1],] <- apply(drugInfo(object)[myx,], 2, paste, collapse="///")
      drugInfo(object) <- drugInfo(object)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(object))>0){
        sensMatch <- match(colnames(sensNumber(object)), rownames(drugInfo(object)))
      }
      if(dim(pertNumber(object))[[2]]>0){
        pertMatch <- match(dimnames(pertNumber(object))[[2]], rownames(drugInfo(object)))
      }
      curMatch <- match(rownames(object@curation$drug),rownames(drugInfo(object)))
    }
  } else {
    if (dim(pertNumber(object))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(object))[[2]], rownames(drugInfo(object)))
    }
    if (ncol(sensNumber(object))>0){
      sensMatch <- match(colnames(sensNumber(object)), rownames(drugInfo(object)))
    }
    curMatch <- match(rownames(object@curation$drug),rownames(drugInfo(object)))
  }

  if (dim(pertNumber(object))[[2]]>0){
    dimnames(pertNumber(object))[[2]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(object))>0){
    colnames(sensNumber(object)) <- new.ids[sensMatch]
  }
  rownames(object@curation$drug) <- new.ids[curMatch]
  rownames(drugInfo(object)) <- new.ids


  return(object)
}

.summarizeSensitivityNumbers <- function(object) {
  
  ## TODO:: Checks don't like assigning to global evnironment. Can we return this?
  assign("object_sumSenNum", object) # Removed envir=.GlobalEnv
  if (object@datasetType != "sensitivity" && object@datasetType != "both") {
    stop ("Data type must be either sensitivity or both")
  }
  
  ## unique drug identifiers
  # drugn <- sort(unique(object@sensitivity$info[ , "drugid"]))
  
  ## consider all drugs
  drugn <- rownames(object@drug)
  
  ## unique drug identifiers
  # celln <- sort(unique(object@sensitivity$info[ , "cellid"]))
  
  ## consider all cell lines
  celln <- rownames(object@cell)
  
  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn), dimnames=list(celln, drugn))
  drugids <- object@sensitivity$info[ , "drugid"]
  cellids <- object@sensitivity$info[ , "cellid"]
  cellids <- cellids[grep("///", drugids, invert=TRUE)]
  drugids <- drugids[grep("///", drugids, invert=TRUE)]
  
  
  tt <- table(cellids, drugids)
  sensitivity.info[rownames(tt), colnames(tt)] <- tt
  
    return(sensitivity.info)
}

#' @importFrom CoreGx .summarizeMolecularNumbers
.summarizeMolecularNumbers <- function(object) {
  CoreGx::.summarizeMolecularNumbers
}

.summarizePerturbationNumbers <- function(object) {

  if (object@datasetType != "perturbation" && object@datasetType != "both") {
    stop ("Data type must be either perturbation or both")
  }
  
  ## unique drug identifiers
  # drugn <- sort(unique(unlist(lapply(object@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "drugid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "drugid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all drugs
  drugn <- rownames(object@drug)
  
  ## unique cell line identifiers
  # celln <- sort(unique(unlist(lapply(object@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "cellid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "cellid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all cell lines
  celln <- rownames(object@cell)
  
  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(object@molecularProfiles)), dimnames=list(celln, drugn, names((object@molecularProfiles))))

    for (i in seq_len(length(object@molecularProfiles))) {
      if (nrow(SummarizedExperiment::colData(object@molecularProfiles[[i]])) > 0 && all(is.element(c("cellid", "drugid"), colnames(SummarizedExperiment::colData(object@molecularProfiles[[i]]))))) {
      tt <- table(SummarizedExperiment::colData(object@molecularProfiles[[i]])[ , "cellid"], SummarizedExperiment::colData(object@molecularProfiles[[i]])[ , "drugid"])
        perturbation.info[rownames(tt), colnames(tt), names(object@molecularProfiles)[i]] <- tt
      }
    }
  
    return(perturbation.info)
}

#' A function to verify the structure of a PharmacoSet
#' 
#' This function checks the structure of a PharamcoSet, ensuring that the
#' correct annotations are in place and all the required slots are filled so
#' that matching of cells and drugs can be properly done across different types
#' of data and with other studies.
#' 
#' @examples
#' data(CCLEsmall)
#' checkPsetStructure(CCLEsmall)
#' 
#' @param object A \code{PharmacoSet} to be verified
#' @param plotDist Should the function also plot the distribution of molecular data?
#' @param result.dir The path to the directory for saving the plots as a string
#' 
#' @return Prints out messages whenever describing the errors found in the 
#'   structure of the object object passed in.
#' 
#' @importFrom graphics hist
#' @importFrom grDevices dev.off pdf
#' 
#' @export
checkPsetStructure <-
  function(object, plotDist=FALSE, result.dir=".") {
    
    # Make directory to store results if it doesn't exist
    if(!file.exists(result.dir) & plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }
    
    #####
    # Checking molecularProfiles
    #####
    # Can this be parallelized or does it mess with the order of printing warnings?
    for( i in seq_along(object@molecularProfiles)) {
      profile <- object@molecularProfiles[[i]]
      nn <- names(object@molecularProfiles)[i]
      
      # Testing plot rendering for rna and rnaseq
      if((S4Vectors::metadata(profile)$annotation == "rna" | S4Vectors::metadata(profile)$annotation == "rnaseq") & plotDist)
      {
        pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
        hist(assays(profile)[[1]], breaks = 100)
        dev.off()
      }
      
      
      ## Test if sample and feature annotations dimensions match the assay
      warning(ifelse(nrow(rowData(profile)) != nrow(assays(profile)[[1]]),
                     sprintf("%s: number of features in fData is different from 
                             SummarizedExperiment slots", nn),
                     sprintf("%s: rowData dimension is OK", nn)
                     )
              )
      warning(ifelse(nrow(colData(profile)) != ncol(assays(profile)[[1]]),
                     sprintf("%s: number of cell lines in pData is different 
                             from expression slots", nn),
                     sprintf("%s: colData dimension is OK", nn)
                     )
              )
      
      
      # Checking sample metadata for required columns
      warning(ifelse("cellid" %in% colnames(colData(profile)), "", 
                     sprintf("%s: cellid does not exist in colData (samples) 
                             columns", nn)))
      warning(ifelse("batchid" %in% colnames(colData(profile)), "", 
                     sprintf("%s: batchid does not exist in colData (samples) 
                             columns", nn)))
      
      # Checking mDataType of the SummarizedExperiment for required columns
      if(S4Vectors::metadata(profile)$annotation == "rna" | 
         S4Vectors::metadata(profile)$annotation == "rnaseq")
      {
        warning(ifelse("BEST" %in% colnames(rowData(profile)), "BEST is OK", 
                       sprintf("%s: BEST does not exist in rowData (features) 
                               columns", nn)))
        warning(ifelse("Symbol" %in% colnames(rowData(profile)), "Symbol is OK", 
                       sprintf("%s: Symbol does not exist in rowData (features) 
                               columns", nn)))
      }
      
      # Check that all cellids from the object are included in molecularProfiles
      if("cellid" %in% colnames(rowData(profile))) {
        if(!all(colData(profile)[,"cellid"] %in% rownames(object@cell))) {
          warning(sprintf("%s: not all the cell lines in this profile are in 
                          cell lines slot", nn))
        }
      }else {
        warning(sprintf("%s: cellid does not exist in colData (samples)", nn))
      }
    }
    
    #####
    # Checking cell
    #####
    if("tissueid" %in% colnames(object@cell)) {
      if("unique.tissueid" %in% colnames(object@curation$tissue))
      {
        if(length(intersect(rownames(object@curation$tissue), 
                            rownames(object@cell))) != nrow(object@cell)) {
          message("rownames of curation tissue slot should be the same as cell 
                  slot (curated cell ids)")
        } else{
          if(length(intersect(object@cell$tissueid, 
                              object@curation$tissue$unique.tissueid)) != 
             length(table(object@cell$tissueid))){
            message("tissueid should be the same as unique tissue id from tissue 
                    curation slot")
          }
        }
      } else {
        message("unique.tissueid which is curated tissue id across data set 
                should be a column of tissue curation slot")
      }
      if(any(is.na(object@cell[,"tissueid"]) | object@cell[,"tissueid"]=="", 
             na.rm=TRUE)){
        message(sprintf("There is no tissue type for this cell line(s): %s", 
                        paste(rownames(object@cell)[which(is.na(
                          object@cell[,"tissueid"]) | 
                            object@cell[,"tissueid"]=="")], collapse=" ")))
      }
    } else {
      warning("tissueid does not exist in cell slot")
    }
    
    if("unique.cellid" %in% colnames(object@curation$cell)) {
      if(length(intersect(object@curation$cell$unique.cellid, 
                          rownames(object@cell))) != nrow(object@cell)) {
        print("rownames of cell slot should be curated cell ids")
      }
    } else {
      print("unique.cellid which is curated cell id across data set should be a 
            column of cell curation slot")
    }
#     if("cellid" %in% colnames(object@cell)) {
#       if(length(intersect(object@curation$cell$cellid, rownames(object@cell)))
#    != nrow(object@cell)) {
#         print("values of cellid column should be curated cell line ids")
#       }
#     } else {
#       print("cellid which is curated cell id across data set should be a column of cell slot")
#     }
    
    if(length(intersect(rownames(object@curation$cell), 
                        rownames(object@cell))) != nrow(object@cell)) {
      print("rownames of curation cell slot should be the same as cell slot 
            (curated cell ids)")
    }
    
    if("unique.drugid" %in% colnames(object@curation$drug)) {
      if(length(intersect(object@curation$drug$unique.drugid, 
                          rownames(object@drug))) != nrow(object@drug)) {
        print("rownames of drug slot should be curated drug ids")
      }
    } else {
      print("unique.drugid which is curated drug id across data set should be a 
            column of drug curation slot")
    }
    
#     if("drugid" %in% colnames(object@drug)) {
#       if(length(intersect(object@curation$drug$drugid, 
#    rownames(object@drug))) != nrow(object@drug)) {
#         print("values of drugid column should be curated drug ids")
#       }
#     } else {
#       print("drugid which is curated drug id across data set should be a 
#    column of drug slot")
#     }
    
    if(length(intersect(rownames(object@curation$cell), 
                        rownames(object@cell))) != nrow(object@cell)) {
      print("rownames of curation drug slot should be the same as drug 
            slot (curated drug ids)")
    }
    
    if(!is(object@cell, "data.frame")) {
      warning("cell slot class type should be dataframe")
    }
    if(!is(object@drug, "data.frame")) {
      warning("drug slot class type should be dataframe")
    }
    if(object@datasetType %in% c("sensitivity", "both"))
    {
      if(!is(object@sensitivity$info, "data.frame")) {
        warning("sensitivity info slot class type should be dataframe")
      }
      if("cellid" %in% colnames(object@sensitivity$info)) {
        if(!all(object@sensitivity$info[,"cellid"] %in% rownames(object@cell))){
          warning("not all the cell lines in sensitivity data are in cell slot")
        }
      }else {
        warning("cellid does not exist in sensitivity info")
      }
      if("drugid" %in% colnames(object@sensitivity$info)) {
        drug.ids <- unique(object@sensitivity$info[,"drugid"])
        drug.ids <- drug.ids[grep("///",drug.ids, invert=TRUE)]
        if(!all(drug.ids %in% rownames(object@drug))) {
          print("not all the drugs in sensitivity data are in drug slot")
        }
      }else {
        warning("drugid does not exist in sensitivity info")
      }
      
      if(any(!is.na(object@sensitivity$raw))) {
        if(!all(dimnames(object@sensitivity$raw)[[1]] %in% 
                rownames(object@sensitivity$info))) {
          warning("For some experiments there is raw sensitivity data but no 
                  experiment information in sensitivity info")
        }
      }
      if(!all(rownames(object@sensitivity$profiles) %in% 
              rownames(object@sensitivity$info))) {
        warning("For some experiments there is sensitivity profiles but no 
                experiment information in sensitivity info")
      }
    }
  }

