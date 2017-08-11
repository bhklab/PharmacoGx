#' A Class to Contain PharmacoGenomic datasets together with their curations
#' 
#' The PharmacoSet (PSet) class was developed to contain and organise large 
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
#' @param pSet A \code{PharmacoSet} object
#' @param mDataType A \code{character} with the type of molecular data to return/update
#' @param object A \code{PharmacoSet} object
#' @param value A replacement value
#' 
#' @slot annotation A \code{list} of annotation data about the PharmacoSet,
#'    including the \code{$name} and the session information for how the object
#'    was creating, detailing the exact versions of R and all the packages used
#' @slot molecularProfiles A \code{list} containing 4 \code{Biobase::ExpressionSet} 
#'   type object for holding data for RNA, DNA, SNP and Copy Number Variation 
#'   measurements respectively, with associated \code{fData} and \code{pData} 
#'   containing the row and column metadata
#' @slot cell A \code{data.frame} containg the annotations for all the cell 
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
#' @return An object of the PharmacoSet class
.PharmacoSet <- setClass("PharmacoSet", slots = list(
                                                     annotation = "list",
                                                     molecularProfiles = "list",
                                                     cell="data.frame", 
                                                     drug="data.frame", 
                                                     datasetType="character", 
                                                     sensitivity="list",
                                                     perturbation="list",
                                                     curation="list"
                                                     # tables="array",
                                                     # table.summary="list",
                                                     # dateCreated="character",
                                                     ))


# The default constructor above does a poor job of explaining the required structure of a PharmacoSet. 
# The constructor function defined below guides the user into providing the required components of the curation and senstivity lists
# and hides the annotation slot which the user does not need to manually fill. 
# This also follows the design of the Expression Set class.

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
#' @param name A \code{character} string detailing the name of the dataset
#' @param molecularProfiles A \code{list} of ExpressionSet objects containing
#'   molecular profiles 
#' @param cell A \code{data.frame} containg the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @param drug A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @param sensitivityInfo A \code{data.frame} containing the information for the
#'   sensitivity experiments
#' @param sensitivityRaw A 3 Dimensional \code{array} contaning the raw drug
#'   dose response data for the sensitivity experiments
#' @param sensitivityProfiles \code{data.frame} containing drug sensitivity profile 
#'   statistics such as IC50 and AUC
#' @param sensitivityN,perturbationN A \code{data.frame} summarizing the
#'   available sensitivity/perturbation data
#' @param curationDrug,curationCell,curationTissue A \code{data.frame} mapping
#'   the names for drugs, cells and tissues used in the data set to universal
#'   identifiers used between different PharmacoSet objects
#' @param datasetType A \code{character} string of 'sensitivity',
#'   'preturbation', or both detailing what type of data can be found in the
#'   PharmacoSet, for proper processing of the data
#' @param verify \code{boolean} Should the function verify the PharmacoSet and
#'   print out any errors it finds after construction?
#' @return An object of class PharmacoSet
#' @export
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
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
                          verify = TRUE)
{
    datasetType <- match.arg(datasetType)
    
    annotation <- list()
    annotation$name <- as.character(name)
    annotation$dateCreated <- date()
    annotation$sessionInfo <- sessionInfo()
    annotation$call <- match.call()
    
    #molecularProfiles <- list("dna"=dna, "rna"=rna, "snp"=snp, "cnv"=cnv)
    for (i in 1:length(molecularProfiles)){
        if (class(molecularProfiles[[i]]) != "ExpressionSet"){
            stop(sprintf("Please provide the %s data as an ExpressionSet", names(molecularProfiles[i])))
        }else{
      Biobase::fData(molecularProfiles[[i]]) <- Biobase::fData(molecularProfiles[[i]])[rownames(Biobase::exprs(molecularProfiles[[i]])), , drop=FALSE]
      Biobase::pData(molecularProfiles[[i]]) <- Biobase::pData(molecularProfiles[[i]])[colnames(Biobase::exprs(molecularProfiles[[i]])), , drop=FALSE]
        }
    
    }
    #if (class(cell)!="data.frame"){
    #    stop("Please provide the cell line annotations as a data frame.")
    #}
    #if (class(drug)!="data.frame"){
    #    stop("Please provide the drug annotations as a data frame.")
    #}
    
    sensitivity <- list()
    
    if (!all(rownames(sensitivityInfo) == rownames(sensitivityProfiles) & rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])){
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
        perturbation$info <- "The metadata for the perturbation experiments is available for each molecular type by calling the appropriate info function. \n For example, for RNA transcriptome perturbations, the metadata can be accessed using rnaInfo(pSet)."
    } else {
        perturbation$info <- "Not a perturbation dataset."
    }
    
    pSet  <- .PharmacoSet(annotation=annotation, molecularProfiles=molecularProfiles, cell=as.data.frame(cell), drug=as.data.frame(drug), datasetType=datasetType, sensitivity=sensitivity, perturbation=perturbation, curation=curation)
    if (verify) { checkPSetStructure(pSet)}
  if(length(sensitivityN) == 0 & datasetType %in% c("sensitivity", "both")) {
    pSet@sensitivity$n <- .summarizeSensitivityNumbers(pSet)
  }
    if(length(perturbationN) == 0  & datasetType %in% c("perturbation", "both")) {
      pSet@perturbation$n <- .summarizePerturbationNumbers(pSet)
    }
  return(pSet)
}
    

#' cellInfo Generic
#' 
#' Generic for cellInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' cellInfo(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} to retrieve cell info from
#' 
#' @return a \code{data.frame} with the cell annotations
setGeneric("cellInfo", function(pSet) standardGeneric("cellInfo"))

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
setGeneric("cellInfo<-", function(object, value) standardGeneric("cellInfo<-"))
#' @describeIn PharmacoSet Returns the annotations for all the cell lines tested on in the PharmacoSet
#' @export
setMethod(cellInfo, "PharmacoSet", function(pSet){
  pSet@cell
})
#' @describeIn PharmacoSet Update the cell line annotations
#' @export
setReplaceMethod("cellInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){
  if(is.null(rownames(value))){
    stop("Please provide the cell_id as rownames for the cell line annotations")
  }
  
  object@cell <- value
  object
})
#' drugInfo Generic
#' 
#' Generic for drugInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' drugInfo(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric("drugInfo", function(pSet) standardGeneric("drugInfo"))

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
#' @return Updated \code{PharmacoSet}
setGeneric("drugInfo<-", function(object, value) standardGeneric("drugInfo<-"))
#' @describeIn PharmacoSet Returns the annotations for all the drugs tested in the PharmacoSet
#' @export
setMethod(drugInfo, "PharmacoSet", function(pSet){

  pSet@drug
})
#' @describeIn PharmacoSet Update the drug annotations
#' @export
setReplaceMethod("drugInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){
  object@drug <- value
  object
})

#' phenoInfo Generic
#' 
#' Generic for phenoInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' phenoInfo(CCLEsmall, mDataType="rna")
#' 
#' @param pSet The \code{PharmacoSet} to retrieve rna annotations from
#' @param mDataType the type of molecular data 
#' @return a \code{data.frame} with the experiment info
setGeneric("phenoInfo", function(pSet, mDataType) standardGeneric("phenoInfo"))
#' @describeIn PharmacoSet Return the experiment info from the given type of molecular data in PharmacoSet 
#' @export
setMethod(phenoInfo, "PharmacoSet", function(pSet, mDataType){
    
  if(mDataType %in% names(pSet@molecularProfiles)){
    return(Biobase::pData(pSet@molecularProfiles[[mDataType]]))}else{
      return(NULL)
    }
    
})

#' phenoInfo<- Generic
#' 
#' Generic for phenoInfo replace method 
#' 
#' @examples
#' 
#' data(CCLEsmall)
#' phenoInfo(CCLEsmall, mDataType="rna") <- phenoInfo(CCLEsmall, mDataType="rna")
#' 
#' @param object The \code{PharmacoSet} to retrieve molecular experiment annotations from
#' @param mDataType the type of molecular data 
#' @param value a \code{data.frame} with the new experiment annotations
#' @return The updated \code{PharmacoSet}
setGeneric("phenoInfo<-", function(object, mDataType, value) standardGeneric("phenoInfo<-"))
#' @describeIn PharmacoSet Update the the given type of molecular data experiment info in the PharmacoSet 
#' @export
setReplaceMethod("phenoInfo", signature = signature(object="PharmacoSet", mDataType ="character",value="data.frame"), function(object, mDataType, value){

  if(mDataType %in% names(object@molecularProfiles)){Biobase::pData(object@molecularProfiles[[mDataType]]) <- value}
    object
})

#' molecularProfiles Generic
#' 
#' Generic for molecularProfiles method 
#' 
#' @examples
#' data(CCLEsmall)
#' molecularProfiles(CCLEsmall, "rna")
#' 
#' @param pSet The \code{PharmacoSet} to retrieve molecular profiles from
#' @param mDataType the type of molecular data 
#' @return a \code{data.frame} with the experiment info
setGeneric("molecularProfiles", function(pSet, mDataType) standardGeneric("molecularProfiles"))
#' @describeIn PharmacoSet Return the given type of molecular data from the PharmacoSet 
#' @export
setMethod(molecularProfiles, "PharmacoSet", function(pSet, mDataType){
    
  if(mDataType %in% names(pSet@molecularProfiles)){
    return(Biobase::exprs(pSet@molecularProfiles[[mDataType]]))}else{
      return(NULL)
    }
    
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
#' @param value A \code{matrix} with the new profiles
#' @return Updated \code{PharmacoSet}
setGeneric("molecularProfiles<-", function(object, mDataType, value) standardGeneric("molecularProfiles<-"))
#' @describeIn PharmacoSet Update the given type of molecular data from the PharmacoSet 
#' @export
setReplaceMethod("molecularProfiles", signature = signature(object="PharmacoSet", mDataType ="character",value="matrix"), function(object, mDataType, value){

  if (mDataType %in% names(object@molecularProfiles)) {
    Biobase::exprs(object@molecularProfiles[[mDataType]]) <- value
  }
    object
})

#' featureInfo Generic
#' 
#' Generic for featureInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' featureInfo(CCLEsmall, "rna")
#' 
#' @param pSet The \code{PharmacoSet} to retrieve feature annotations from
#' @param mDataType the type of molecular data 
#' @return a \code{data.frame} with the experiment info
setGeneric("featureInfo", function(pSet, mDataType) standardGeneric("featureInfo"))
#' @describeIn PharmacoSet Return the feature info for the given molecular data 
#' @export
setMethod(featureInfo, "PharmacoSet", function(pSet, mDataType){
  if(mDataType %in% names(pSet@molecularProfiles)){
    return(Biobase::fData(pSet@molecularProfiles[[mDataType]]))}else{
      return(NULL)
    }
  
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
#' @param value A \code{data.frame} with the new feature annotations
#' @return Updated \code{PharmacoSet}
setGeneric("featureInfo<-", function(object, mDataType, value) standardGeneric("featureInfo<-"))
#' @describeIn PharmacoSet Replace the gene info for the molecular data
#' @export
setReplaceMethod("featureInfo", signature = signature(object="PharmacoSet", mDataType ="character",value="data.frame"), function(object, mDataType, value){
  
  if(mDataType %in% names(object@molecularProfiles)){Biobase::fData(object@molecularProfiles[[mDataType]]) <- value}
  
  object
})

#' sensitivityInfo Generic
#' 
#' Generic for sensitivityInfo method 
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} to retrieve sensitivity experiment annotations from
#' @return a \code{data.frame} with the experiment info
setGeneric("sensitivityInfo", function(pSet) standardGeneric("sensitivityInfo"))
#' @describeIn PharmacoSet Return the drug dose sensitivity experiment info
#' @export
setMethod(sensitivityInfo, "PharmacoSet", function(pSet){
    
    return(pSet@sensitivity$info)
    
})

#' sensitivityInfo<- Generic
#' 
#' A generic for the sensitivityInfo replacement method
#' 
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityInfo(CCLEsmall) <- sensitivityInfo(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{data.frame} with the new sensitivity annotations
#' @return Updated \code{PharmacoSet} 
setGeneric("sensitivityInfo<-", function(object, value) standardGeneric("sensitivityInfo<-"))
#' @describeIn PharmacoSet Update the sensitivity experiment info
#' @export
setReplaceMethod("sensitivityInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

    object@sensitivity$info <- value
    object
})


#' sensitivityProfiles Generic
#' 
#' Generic for sensitivityProfiles method 
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityProfiles(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} to retrieve sensitivity experiment data from
#' @return a \code{data.frame} with the experiment info
setGeneric("sensitivityProfiles", function(pSet) standardGeneric("sensitivityProfiles"))
#' @describeIn PharmacoSet Return the phenotypic data for the drug dose sensitivity
#' @export
setMethod(sensitivityProfiles, "PharmacoSet", function(pSet){
    
    return(pSet@sensitivity$profiles)
    
})

#' sensitivityProfiles<- Generic
#' 
#' A generic for the sensitivityProfiles replacement method
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityProfiles(CCLEsmall) <- sensitivityProfiles(CCLEsmall)
#' 
#' @param object The \code{PharmacoSet} to update
#' @param value A \code{data.frame} with the new sensitivity profiles. If a matrix object is passed in, converted to data.frame before assignment
#' @return Updated \code{PharmacoSet} 
setGeneric("sensitivityProfiles<-", function(object, value) standardGeneric("sensitivityProfiles<-"))
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

    object@sensitivity$profiles <- value
    object
})
#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityProfiles", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

    object@sensitivity$profiles <- as.data.frame(value)
    object
})
#' sensitivityMeasures Generic
#' 
#' A generic for the sensitivityMeasures  method
#' 
#' @examples
#' data(CCLEsmall)
#' sensitivityMeasures(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} 
#' @return A \code{character} vector of all the available sensitivity measures
setGeneric("sensitivityMeasures", function(pSet) standardGeneric("sensitivityMeasures"))
#' @describeIn PharmacoSet Returns the available sensitivity profile
#'   summaries, for example, whether there are IC50 values available
#' @export
setMethod(sensitivityMeasures, "PharmacoSet", function(pSet){
    
    return(colnames(sensitivityProfiles(pSet)))
    
})

#' drugNames Generic
#' 
#' A generic for the drugNames method
#' 
#' @examples
#' data(CCLEsmall)
#' drugNames(CCLEsmall)
#' 
#' @param pSet The \code{PharmacoSet} to return drug names from
#' @return A vector of the drug names used in the PharmacoSet
setGeneric("drugNames", function(pSet) standardGeneric("drugNames"))
#' @describeIn PharmacoSet Return the names of the drugs used in the PharmacoSet
#' @export
setMethod(drugNames, "PharmacoSet", function(pSet){
  
  # if (unique){
#     unique(pData(pSet)[["drugid"]])
#   } else {
#     pData(pSet)[["drugid"]]
#   
#  }
  rownames(drugInfo(pSet))

})

#' drugNames<- Generic
#' 
#' A generic for the drugNames replacement method
#' 
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
#' @param pSet The \code{PharmacoSet} to return cell names from
#' @return A vector of the cell names used in the PharmacoSet
setGeneric("cellNames", function(pSet) standardGeneric("cellNames"))
#' @describeIn PharmacoSet Return the cell names used in the dataset
#' @export
setMethod(cellNames, "PharmacoSet", function(pSet){
  
  rownames(cellInfo(pSet))
  
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
#' @return Updated \code{PharmacoSet} 
setGeneric("cellNames<-", function(object, value) standardGeneric("cellNames<-"))
#' @describeIn PharmacoSet Update the cell names used in the dataset
#' @export
setReplaceMethod("cellNames", signature = signature(object="PharmacoSet",value="character"), function(object, value){
    
    object <- updateCellId(object, value)
    return(object)
    })


    #### TODO:: set replace method for genenames
#' fNames Generic
#' 
#' A generic for the fNames method
#' 
#' @examples
#' data(CCLEsmall)
#' fNames(CCLEsmall, "rna")
#' 
#' @param pSet The \code{PharmacoSet} 
#' @param mDataType The molecular data type to return feature names for
#' @return A \code{character} vector of the feature names
setGeneric("fNames", function(pSet, mDataType) standardGeneric("fNames"))
#' @describeIn PharmacoSet Return the feature names used in the dataset
#' @export
setMethod(fNames, "PharmacoSet", function(pSet, mDataType){
  if (mDataType %in% names(pSet@molecularProfiles)) {
    rownames(featureInfo(pSet, mDataType))
  } else {
    stop("Molecular data type name specified is not part of this PharmacoSet")
  }
})

#' dateCreated Generic
#' 
#' A generic for the dateCreated method
#' 
#' @examples
#' data(CCLEsmall)
#' dateCreated(CCLEsmall)
#' 
#' @param pSet A \code{PharmacoSet} 
#' @return The date the PharmacoSet was created
setGeneric("dateCreated", function(pSet) standardGeneric("dateCreated"))
#' @describeIn PharmacoSet Return the date the PharmacoSet was created
#' @export
setMethod(dateCreated, "PharmacoSet", function(pSet) {
  pSet@annotation$dateCreated
})

#' pSetName Generic
#' 
#' A generic for the pSetName method
#' 
#' @examples
#' data(CCLEsmall)
#' pSetName(CCLEsmall)
#' 
#' @param pSet A \code{PharmacoSet} 
#' @return The name of the PharmacoSet
setGeneric("pSetName", function(pSet) standardGeneric("pSetName"))
#' @describeIn PharmacoSet Return the name of the PharmacoSet 
#' @export
setMethod(pSetName, "PharmacoSet", function(pSet){
    
    return(pSet@annotation$name)
    
})

#' pertNumber Generic
#' 
#' A generic for the pertNumber method
#' 
#' @examples
#' data(CCLEsmall)
#' pertNumber(CCLEsmall)
#' 
#' @param pSet A \code{PharmacoSet} 
#' @return A 3D \code{array} with the number of perturbation experiments per drug and cell line, and data type
setGeneric("pertNumber", function(pSet) standardGeneric("pertNumber"))
#' @describeIn PharmacoSet Return the summary of available perturbation
#'   experiments
#' @export
setMethod(pertNumber, "PharmacoSet", function(pSet){
    
    return(pSet@perturbation$n)
    
})


#' sensNumber Generic
#' 
#' A generic for the sensNumber method
#' 
#' @examples
#' data(CCLEsmall)
#' sensNumber(CCLEsmall)
#' 
#' @param pSet A \code{PharmacoSet} 
#' @return A \code{data.frame} with the number of sensitivity experiments per drug and cell line
setGeneric("sensNumber", function(pSet) standardGeneric("sensNumber"))
#' @describeIn PharmacoSet Return the summary of available sensitivity
#'   experiments
#' @export
setMethod(sensNumber, "PharmacoSet", function(pSet){
  
  return(pSet@sensitivity$n)
  
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
#' @return The updated \code{PharmacoSet} 
setGeneric("pertNumber<-", function(object, value) standardGeneric("pertNumber<-"))
#' @describeIn PharmacoSet Update the summary of available perturbation
#'   experiments
#' @export
setReplaceMethod('pertNumber', signature = signature(object="PharmacoSet",value="array"), function(object, value){
  
  object@perturbation$n <- value
  object
  
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
#' @return The updated \code{PharmacoSet} 
setGeneric("sensNumber<-", function(object, value) standardGeneric("sensNumber<-"))
#' @describeIn PharmacoSet Update the summary of available sensitivity
#'   experiments
#' @export
setReplaceMethod('sensNumber', signature = signature(object="PharmacoSet",value="matrix"), function(object, value){
  
  object@sensitivity$n <- value
  object
  
})

#' Show a PharamcoSet
#' 
#' @param object \code{PharmacoSet}
#' 
#' @examples
#' data(CCLEsmall)
#' CCLEsmall
#' 
#' @return Prints the PharmacoSet object to the output stream, and returns invisible NULL. 
#' @export
setMethod("show", signature=signature(object="PharmacoSet"), 
    function(object) {
        cat("Name: ", pSetName(object), "\n")
        cat("Date Created: ", dateCreated(object), "\n")
    cat("Number of cell lines: ", nrow(cellInfo(object)), "\n")
    cat("Number of drug compounds: ", nrow(drugInfo(object)), "\n")
        if("dna" %in% names(object@molecularProfiles)){cat("DNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="dna")), "\n")}
      if("rna" %in% names(object@molecularProfiles)){cat("RNA: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rna")), "\n")}
      if("rnaseq" %in% names(object@molecularProfiles)){cat("RNASeq: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="rnaseq")), "\n")}
      if("snp" %in% names(object@molecularProfiles)){cat("SNP: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="snp")), "\n")}
      if("cnv" %in% names(object@molecularProfiles)){cat("CNV: \n");cat("\tDim: ", dim(molecularProfiles(object, mDataType="cnv")), "\n")}
        cat("Drug pertubation: \n")
        cat("\tPlease look at pertNumber(pSet) to determine number of experiments for each drug-cell combination.\n")
        cat("Drug sensitivity: \n")
        cat("\tNumber of Experiments: ",nrow(sensitivityInfo(object)),"\n")
        cat("\tPlease look at sensNumber(pSet) to determine number of experiments for each drug-cell combination.\n")
    })


#' mDataNames
#' 
#' Returns the molecular data names for the PharmacoSet.
#' 
#' @examples
#' data(CCLEsmall)
#' mDataNames(CCLEsmall)
#' 
#' @param pSet PharamcoSet object
#' @return Vector of names of the molecular data types
#' @export
mDataNames <- function(pSet){

  return(names(pSet@molecularProfiles))

}

#'`[`
#'
#'@param x PSet
#'@param i Cell lines to keep in PSet
#'@param j Drugs to keep in PSet
#'@param ... further arguments
#'@param drop A boolean flag of whether to drop single dimensions or not
#'@return Returns the subsetted PSet
#'@export
setMethod(`[`, "PharmacoSet", function(x, i, j, ..., drop = FALSE){
  if(is.character(i)&&is.character(j)){
    return(subsetTo(x, cells=i, drugs=j,  molecular.data.cells=i))
  } 
  else if(is.numeric(i) && is.numeric(j) && (as.integer(i)==i) && (as.integer(j)==j)){
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
#' PSet <- subsetTo(CCLEsmall,drugs = CCLEdrugs[1], cells = CCLEcells[1])
#' PSet
#' 
#' @param pSet A \code{PharmacoSet} to be subsetted
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
#' @return A PharmacoSet with only the selected drugs and cells
#' @export
# subsetTo <- function(pSet, cells=NULL, drugs=NULL, exps=NULL, molecular.data.cells=NULL, keep.controls=TRUE) {
subsetTo <- function(pSet, cells=NULL, drugs=NULL, molecular.data.cells=NULL, keep.controls=TRUE, ...) {
  drop=FALSE
  
  adArgs = list(...)
  if ("exps" %in% names(adArgs)) {
  	exps <- adArgs[["exps"]]
  	if(class(exps)=="data.frame"){
  		exps2 <- exps[[pSetName(pSet)]]
  		names(exps2) <- rownames(exps)
  		exps <- exps2
  	} else{
  		exps <- exps[[pSetName(pSet)]]
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
  
    ### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
    ### function in lapply, so it does not recognize them as missing
  
  pSet@molecularProfiles <- lapply(pSet@molecularProfiles, function(eset, cells, drugs, molecular.data.cells){
    
    molecular.data.type <- ifelse(length(grep("rna", Biobase::annotation(eset)) > 0), "rna", Biobase::annotation(eset))
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cells <- molecular.data.cells[[molecular.data.type]]
    }
      column_indices <- NULL
  
      if (length(cells)==0 && length(drugs)==0) {
          column_indices <- 0:ncol(eset)
      }
      if(length(cells)==0 && pSet@datasetType=="sensitivity") {
        column_indices <- 0:ncol(eset)
      }
  
      cell_line_index <- NULL
      if(length(cells)!=0) {
        if (!all(cells %in% cellNames(pSet))) {
              stop("Some of the cell names passed to function did not match to names in the PharmacoSet. Please ensure you are using cell names as returned by the cellNames function")
        }
          cell_line_index <- which(Biobase::pData(eset)[["cellid"]] %in% cells)
        # if (length(na.omit(cell_line_index))==0){
    #       stop("No cell lines matched")
    #     }
      }
      drugs_index <- NULL
      if(pSet@datasetType=="perturbation" || pSet@datasetType=="both"){
        if(length(drugs) != 0) {
            if (!all(drugs %in% drugNames(pSet))){
                  stop("Some of the drug names passed to function did not match to names in the PharmacoSet. Please ensure you are using drug names as returned by the drugNames function")
            }
          drugs_index <- which(Biobase::pData(eset)[["drugid"]] %in% drugs)
          # if (length(drugs_index)==0){
    #         stop("No drugs matched")
    #       }
          if(keep.controls) {
            control_indices <- which(Biobase::pData(eset)[["xptype"]]=="control")
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
  
      row_indices <- 0:nrow(Biobase::exprs(eset))
  
      eset <- eset[row_indices,column_indices]
      return(eset)

  }, cells=cells, drugs=drugs, molecular.data.cells=molecular.data.cells)  
  
  if ((pSet@datasetType == "sensitivity" | pSet@datasetType == "both") & length(exps) != 0) {
      pSet@sensitivity$info <- pSet@sensitivity$info[exps, , drop=drop]
      rownames(pSet@sensitivity$info) <- names(exps)
      if(length(pSet@sensitivity$raw) > 0) {
        pSet@sensitivity$raw <- pSet@sensitivity$raw[exps, , , drop=drop]
        dimnames(pSet@sensitivity$raw)[[1]] <- names(exps)
      }
      pSet@sensitivity$profiles <- pSet@sensitivity$profiles[exps, , drop=drop]
      rownames(pSet@sensitivity$profiles) <- names(exps)
      
      pSet@sensitivity$n <- .summarizeSensitivityNumbers(pSet)
  }
  else if ((pSet@datasetType == "sensitivity" | pSet@datasetType == "both") & (length(drugs) != 0 | length(cells) != 0)) {
    
        drugs_index <- which (sensitivityInfo(pSet)[, "drugid"] %in% drugs)
        cell_line_index <- which (sensitivityInfo(pSet)[,"cellid"] %in% cells)
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
        pSet@sensitivity[names(pSet@sensitivity)[names(pSet@sensitivity)!="n"]] <- lapply(pSet@sensitivity[names(pSet@sensitivity)[names(pSet@sensitivity)!="n"]], function(x,i, drop){
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
		if(pSet@datasetType == "sensitivity" | pSet@datasetType == "both"){
			drugs <- unique(sensitivityInfo(pSet)[["drugid"]])
		}
		if(pSet@datasetType == "perturbation" | pSet@datasetType == "both"){
			drugs <- union(drugs, na.omit(unionList(lapply(pSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["drugid"]])}))))
		}
	}
	if (length(cells)==0) {
		cells <- union(cells, na.omit(unionList(lapply(pSet@molecularProfiles, function(eSet){unique(Biobase::pData(eSet)[["cellid"]])}))))
        if (pSet@datasetType =="sensitivity" | pSet@datasetType == "both"){
            cells <- union(cells, sensitivityInfo(pSet)[["cellid"]])
        }
	}
	drugInfo(pSet) <- drugInfo(pSet)[drugs , , drop=drop]
	cellInfo(pSet) <- cellInfo(pSet)[cells , , drop=drop]
	pSet@curation$drug <- pSet@curation$drug[drugs , , drop=drop]
	pSet@curation$cell <- pSet@curation$cell[cells , , drop=drop]
	pSet@curation$tissue <- pSet@curation$tissue[cells , , drop=drop]
	if (pSet@datasetType == "sensitivity" | pSet@datasetType == "both"  & length(exps) == 0) {
	  pSet@sensitivity$n <- pSet@sensitivity$n[cells, drugs , drop=drop]
	}
	if (pSet@datasetType == "perturbation" | pSet@datasetType == "both") {
	    pSet@perturbation$n <- pSet@perturbation$n[cells,drugs, , drop=drop]
    }
      return(pSet)
}
### TODO:: Add updating of sensitivity Number tables
updateCellId <- function(pSet, new.ids = vector("character")){
  
  if (length(new.ids)!=nrow(cellInfo(pSet))){
    stop("Wrong number of cell identifiers")
  }

  if(pSet@datasetType=="sensitivity"|pSet@datasetType=="both"){
    myx <- match(sensitivityInfo(pSet)[,"cellid"],rownames(cellInfo(pSet)))
    sensitivityInfo(pSet)[,"cellid"] <- new.ids[myx]

  }
  
  
  pSet@molecularProfiles <- lapply(pSet@molecularProfiles, function(eset){
          
      myx <- match(Biobase::pData(eset)[["cellid"]],rownames(cellInfo(pSet)))
      Biobase::pData(eset)[["cellid"]]  <- new.ids[myx]
      return(eset)
        })





  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateCellId. Merging old ids into the same identifier")
    
    if(ncol(sensNumber(pSet))>0){
      sensMatch <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
    }
    if(dim(pertNumber(pSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
    }
    curMatch <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(pSet))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(pSet)[myx[1],] <- apply(sensNumber(pSet)[myx,], 2, sum)
        sensNumber(pSet) <- sensNumber(pSet)[-myx[-1],]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(pSet))[[1]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(pSet)[myx[1],,] <- apply(pertNumber(pSet)[myx,,], c(1,3), sum)
        pertNumber(pSet) <- pertNumber(pSet)[-myx[-1],,]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      pSet@curation$cell[myx[1],] <- apply(pSet@curation$cell[myx,], 2, paste, collapse="///")
      pSet@curation$cell <- pSet@curation$cell[-myx[-1],]
      pSet@curation$tissue[myx[1],] <- apply(pSet@curation$tissue[myx,], 2, paste, collapse="///")
      pSet@curation$tissue <- pSet@curation$tissue[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      cellInfo(pSet)[myx[1],] <- apply(cellInfo(pSet)[myx,], 2, paste, collapse="///")
      cellInfo(pSet) <- cellInfo(pSet)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(pSet))>0){
        sensMatch <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
      }
      if(dim(pertNumber(pSet))[[1]]>0){
        pertMatch <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
      }
      curMatch <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))
    }
  } else {
    if (dim(pertNumber(pSet))[[1]]>0){
      pertMatch <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
    }
    if (ncol(sensNumber(pSet))>0){
      sensMatch <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
    }
    curMatch <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))
  }

  if (dim(pertNumber(pSet))[[1]]>0){
    dimnames(pertNumber(pSet))[[1]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(pSet))>0){
    rownames(sensNumber(pSet)) <- new.ids[sensMatch]
  }
  rownames(pSet@curation$cell) <- new.ids[curMatch]
  rownames(pSet@curation$tissue) <- new.ids[curMatch]
  rownames(cellInfo(pSet)) <- new.ids





  # myx <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))
  # rownames(pSet@curation$cell) <- new.ids[myx]
  # rownames(pSet@curation$tissue) <- new.ids[myx]
  # if (dim(pertNumber(pSet))[[1]]>0){
  #   myx <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
  #   dimnames(pertNumber(pSet))[[1]] <- new.ids[myx]
  # }
  # if (nrow(sensNumber(pSet))>0){
  #   myx <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
  #   rownames(sensNumber(pSet)) <- new.ids[myx]
  # }
  # rownames(cellInfo(pSet)) <- new.ids
  return(pSet)

}

# updateFeatureNames <- function(pSet, new.ids = vector("character")){
#
#   if (length(new.ids)!=nrow(cellInfo(pSet))){
#     stop("Wrong number of cell identifiers")
#   }
#
#   if(pSet@datasetType=="sensitivity"|pSet@datasetType=="both"){
#     myx <- match(sensitivityInfo(pSet)[,"cellid"],rownames(cellInfo(pSet)))
#     sensitivityInfo(pSet)[,"cellid"] <- new.ids[myx]
#
#   }
#
#   pSet@molecularProfiles <- lapply(pSet@molecularProfiles, function(eset){
#
#     myx <- match(pData(eset)[["cellid"]],rownames(cellInfo(pSet)))
#     pData(eset)[["cellid"]]  <- new.ids[myx]
#     return(eset)
#       })
#   myx <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))
#   rownames(pSet@curation$cell) <- new.ids[myx]
#   rownames(pSet@curation$tissue) <- new.ids[myx]
#   if (dim(pertNumber(pSet))[[1]]>0){
#     myx <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
#     dimnames(pertNumber(pSet))[[1]] <- new.ids[myx]
#   }
#   if (nrow(sensNumber(pSet))>0){
#     myx <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
#     rownames(sensNumber(pSet)) <- new.ids[myx]
#   }
#   rownames(cellInfo(pSet)) <- new.ids
#   return(pSet)
#
# }

### TODO:: Add updating of sensitivity Number tables
updateDrugId <- function(pSet, new.ids = vector("character")){

  if (length(new.ids)!=nrow(drugInfo(pSet))){
    stop("Wrong number of drug identifiers")
  }

  if(pSet@datasetType=="sensitivity"|pSet@datasetType=="both"){
    myx <- match(sensitivityInfo(pSet)[,"drugid"],rownames(drugInfo(pSet)))
    sensitivityInfo(pSet)[,"drugid"] <- new.ids[myx]

  }
  if(pSet@datasetType=="perturbation"|pSet@datasetType=="both"){
    pSet@molecularProfiles <- lapply(pSet@molecularProfiles, function(eset){

      myx <- match(Biobase::pData(eset)[["drugid"]],rownames(drugInfo(pSet)))
      Biobase::pData(eset)[["drugid"]]  <- new.ids[myx]
      return(eset)
    })
  }
  

  if(any(duplicated(new.ids))){
    warning("Duplicated ids passed to updateDrugId. Merging old ids into the same identifier")
    
    if(ncol(sensNumber(pSet))>0){
      sensMatch <- match(colnames(sensNumber(pSet)), rownames(drugInfo(pSet)))
    }
    if(dim(pertNumber(pSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(pSet))[[2]], rownames(drugInfo(pSet)))
    }
    curMatch <- match(rownames(pSet@curation$drug),rownames(drugInfo(pSet)))

    duplId <- unique(new.ids[duplicated(new.ids)])
    for(id in duplId){

      if (ncol(sensNumber(pSet))>0){
        myx <- which(new.ids[sensMatch] == id)
        sensNumber(pSet)[,myx[1]] <- apply(sensNumber(pSet)[,myx], 1, sum)
        sensNumber(pSet) <- sensNumber(pSet)[,-myx[-1]]
        # sensMatch <- sensMatch[-myx[-1]]
      }
      if (dim(pertNumber(pSet))[[2]]>0){
        myx <- which(new.ids[pertMatch] == id)
        pertNumber(pSet)[,myx[1],] <- apply(pertNumber(pSet)[,myx,], c(1,3), sum)
        pertNumber(pSet) <- pertNumber(pSet)[,-myx[-1],]
        # pertMatch <- pertMatch[-myx[-1]]
      }

      myx <- which(new.ids[curMatch] == id)
      pSet@curation$drug[myx[1],] <- apply(pSet@curation$drug[myx,], 2, paste, collapse="///")
      pSet@curation$drug <- pSet@curation$drug[-myx[-1],]
      # curMatch <- curMatch[-myx[-1]]

      myx <- which(new.ids == id)
      drugInfo(pSet)[myx[1],] <- apply(drugInfo(pSet)[myx,], 2, paste, collapse="///")
      drugInfo(pSet) <- drugInfo(pSet)[-myx[-1],]
      new.ids <- new.ids[-myx[-1]]
      if(ncol(sensNumber(pSet))>0){
        sensMatch <- match(colnames(sensNumber(pSet)), rownames(drugInfo(pSet)))
      }
      if(dim(pertNumber(pSet))[[2]]>0){
        pertMatch <- match(dimnames(pertNumber(pSet))[[2]], rownames(drugInfo(pSet)))
      }
      curMatch <- match(rownames(pSet@curation$drug),rownames(drugInfo(pSet)))
    }
  } else {
    if (dim(pertNumber(pSet))[[2]]>0){
      pertMatch <- match(dimnames(pertNumber(pSet))[[2]], rownames(drugInfo(pSet)))
    }
    if (ncol(sensNumber(pSet))>0){
      sensMatch <- match(colnames(sensNumber(pSet)), rownames(drugInfo(pSet)))
    }
    curMatch <- match(rownames(pSet@curation$drug),rownames(drugInfo(pSet)))
  }

  if (dim(pertNumber(pSet))[[2]]>0){
    dimnames(pertNumber(pSet))[[2]] <- new.ids[pertMatch]
  }
  if (ncol(sensNumber(pSet))>0){
    colnames(sensNumber(pSet)) <- new.ids[sensMatch]
  }
  rownames(pSet@curation$drug) <- new.ids[curMatch]
  rownames(drugInfo(pSet)) <- new.ids


  return(pSet)
}

.summarizeSensitivityNumbers <- function(pSet) {

  if (pSet@datasetType != "sensitivity" && pSet@datasetType != "both") {
    stop ("Data type must be either sensitivity or both")
  }
  
  ## unique drug identifiers
  # drugn <- sort(unique(pSet@sensitivity$info[ , "drugid"]))
  
  ## consider all drugs
  drugn <- rownames(pSet@drug)
  
  ## unique drug identifiers
  # celln <- sort(unique(pSet@sensitivity$info[ , "cellid"]))
  
  ## consider all cell lines
  celln <- rownames(pSet@cell)
  
  sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn), dimnames=list(celln, drugn))
  drugids <- pSet@sensitivity$info[ , "drugid"]
  cellids <- pSet@sensitivity$info[ , "cellid"]
  cellids <- cellids[grep("///", drugids, invert=TRUE)]
  drugids <- drugids[grep("///", drugids, invert=TRUE)]
  
  
  tt <- table(cellids, drugids)
  sensitivity.info[rownames(tt), colnames(tt)] <- tt
  
    return(sensitivity.info)
}


.summarizeMolecularNumbers <- function(pSet) {
  
  ## consider all molecular types
  mDT <- mDataNames(pSet)
  
  ## consider all cell lines
  celln <- rownames(pSet@cell)
  
  molecular.info <- matrix(0, nrow=length(celln), ncol=length(mDT), dimnames=list(celln, mDT))
  
  for(mDataType in mDT) {
    tt <- table(phenoInfo(pSet, mDataType)$cellid)
    molecular.info[names(tt), mDataType] <- tt

  }
  return(molecular.info)
}


.summarizePerturbationNumbers <- function(pSet) {

  if (pSet@datasetType != "perturbation" && pSet@datasetType != "both") {
    stop ("Data type must be either perturbation or both")
  }
  
  ## unique drug identifiers
  # drugn <- sort(unique(unlist(lapply(pSet@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "drugid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "drugid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all drugs
  drugn <- rownames(pSet@drug)
  
  ## unique cell line identifiers
  # celln <- sort(unique(unlist(lapply(pSet@molecularProfiles, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "cellid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "cellid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all cell lines
  celln <- rownames(pSet@cell)
  
  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(pSet@molecularProfiles)), dimnames=list(celln, drugn, names((pSet@molecularProfiles))))

    for (i in 1:length(pSet@molecularProfiles)) {
      if (nrow(Biobase::pData(pSet@molecularProfiles[[i]])) > 0 && all(is.element(c("cellid", "drugid"), colnames(Biobase::pData(pSet@molecularProfiles[[i]]))))) {
      tt <- table(Biobase::pData(pSet@molecularProfiles[[i]])[ , "cellid"], Biobase::pData(pSet@molecularProfiles[[i]])[ , "drugid"])
        perturbation.info[rownames(tt), colnames(tt), names(pSet@molecularProfiles)[i]] <- tt
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
#' 
#' checkPSetStructure(CCLEsmall)
#' 
#' @param pSet A \code{PharmacoSet} to be verified
#' @param plotDist Should the function also plot the distribution of molecular data?
#' @param result.dir The path to the directory for saving the plots as a string
#' @return Prints out messages whenever describing the errors found in the structure of the pset object passed in. 
#' @export
#' @importFrom graphics hist
#' @importFrom grDevices dev.off pdf

checkPSetStructure <-
  function(pSet, plotDist=FALSE, result.dir=".") {
    if(!file.exists(result.dir) & plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }
    for( i in 1:length(pSet@molecularProfiles)) {
      profile <- pSet@molecularProfiles[[i]]
      nn <- names(pSet@molecularProfiles)[i]
      if((Biobase::annotation(profile) == "rna" | Biobase::annotation(profile) == "rnaseq") & plotDist)
      {
        pdf(file=file.path(result.dir, sprintf("%s.pdf", nn)))
        hist(Biobase::exprs(profile), breaks = 100)
        dev.off()
      }
      warning(ifelse(nrow(Biobase::fData(profile)) != nrow(Biobase::exprs(profile)), sprintf("%s: number of features in fData is different from expression slots", nn), sprintf("%s: fData dimension is OK", nn)))
      warning(ifelse(nrow(Biobase::pData(profile)) != ncol(Biobase::exprs(profile)), sprintf("%s: number of cell lines in pData is different from expression slots", nn), sprintf("%s: pData dimension is OK", nn)))
      warning(ifelse("cellid" %in% colnames(Biobase::pData(profile)), "", sprintf("%s: cellid does not exist in pData columns", nn)))
      warning(ifelse("batchid" %in% colnames(Biobase::pData(profile)), "", sprintf("%s: batchid does not exist in pData columns", nn)))
      if(Biobase::annotation(profile) == "rna" | Biobase::annotation(profile) == "rnaseq")
      {
        warning(ifelse("BEST" %in% colnames(Biobase::fData(profile)), "BEST is OK", sprintf("%s: BEST does not exist in fData columns", nn)))
        warning(ifelse("Symbol" %in% colnames(Biobase::fData(profile)), "Symbol is OK", sprintf("%s: Symbol does not exist in fData columns", nn)))
      }
      if("cellid" %in% colnames(Biobase::pData(profile))) {
        if(!all(Biobase::pData(profile)[,"cellid"] %in% rownames(pSet@cell))) {
          warning(sprintf("%s: not all the cell lines in this profile are in cell lines slot", nn))
        }
      }else {
        warning(sprintf("%s: cellid does not exist in pData", nn))
      }
    }
    if("tissueid" %in% colnames(pSet@cell)) {
      if("unique.tissueid" %in% colnames(pSet@curation$tissue))
      {
        if(length(intersect(rownames(pSet@curation$tissue), rownames(pSet@cell))) != nrow(pSet@cell)) {
          message("rownames of curation tissue slot should be the same as cell slot (curated cell ids)")
        } else{
          if(length(intersect(pSet@cell$tissueid, pSet@curation$tissue$unique.tissueid)) != length(table(pSet@cell$tissueid))){
            message("tissueid should be the same as unique tissue id from tissue curation slot")
          }
        }
      } else {
        message("unique.tissueid which is curated tissue id across data set should be a column of tissue curation slot")
      }
      if(any(is.na(pSet@cell[,"tissueid"]) | pSet@cell[,"tissueid"]=="", na.rm=TRUE)){
        message(sprintf("There is no tissue type for this cell line(s): %s", paste(rownames(pSet@cell)[which(is.na(pSet@cell[,"tissueid"]) | pSet@cell[,"tissueid"]=="")], collapse=" ")))
      }
    } else {
      warning("tissueid does not exist in cell slot")
    }
    
    if("unique.cellid" %in% colnames(pSet@curation$cell)) {
      if(length(intersect(pSet@curation$cell$unique.cellid, rownames(pSet@cell))) != nrow(pSet@cell)) {
        print("rownames of cell slot should be curated cell ids")
      }
    } else {
      print("unique.cellid which is curated cell id across data set should be a column of cell curation slot")
    }
#     if("cellid" %in% colnames(pSet@cell)) {
#       if(length(intersect(pSet@curation$cell$cellid, rownames(pSet@cell))) != nrow(pSet@cell)) {
#         print("values of cellid column should be curated cell line ids")
#       }
#     } else {
#       print("cellid which is curated cell id across data set should be a column of cell slot")
#     }
    
    if(length(intersect(rownames(pSet@curation$cell), rownames(pSet@cell))) != nrow(pSet@cell)) {
      print("rownames of curation cell slot should be the same as cell slot (curated cell ids)")
    }
    
    if("unique.drugid" %in% colnames(pSet@curation$drug)) {
      if(length(intersect(pSet@curation$drug$unique.drugid, rownames(pSet@drug))) != nrow(pSet@drug)) {
        print("rownames of drug slot should be curated drug ids")
      }
    } else {
      print("unique.drugid which is curated drug id across data set should be a column of drug curation slot")
    }
    
#     if("drugid" %in% colnames(pSet@drug)) {
#       if(length(intersect(pSet@curation$drug$drugid, rownames(pSet@drug))) != nrow(pSet@drug)) {
#         print("values of drugid column should be curated drug ids")
#       }
#     } else {
#       print("drugid which is curated drug id across data set should be a column of drug slot")
#     }
    
    if(length(intersect(rownames(pSet@curation$cell), rownames(pSet@cell))) != nrow(pSet@cell)) {
      print("rownames of curation drug slot should be the same as drug slot (curated drug ids)")
    }
    
    if(class(pSet@cell) != "data.frame") {
      warning("cell slot class type should be dataframe")
    }
    if(class(pSet@drug) != "data.frame") {
      warning("drug slot class type should be dataframe")
    }
    if(pSet@datasetType %in% c("sensitivity", "both"))
    {
      if(class(pSet@sensitivity$info) != "data.frame") {
        warning("sensitivity info slot class type should be dataframe")
      }
      if("cellid" %in% colnames(pSet@sensitivity$info)) {
        if(!all(pSet@sensitivity$info[,"cellid"] %in% rownames(pSet@cell))) {
          warning("not all the cell lines in sensitivity data are in cell slot")
        }
      }else {
        warning("cellid does not exist in sensitivity info")
      }
      if("drugid" %in% colnames(pSet@sensitivity$info)) {
        drug.ids <- unique(pSet@sensitivity$info[,"drugid"])
        drug.ids <- drug.ids[grep("///",drug.ids, invert=TRUE)]
        if(!all(drug.ids %in% rownames(pSet@drug))) {
          print("not all the drugs in sensitivity data are in drug slot")
        }
      }else {
        warning("drugid does not exist in sensitivity info")
      }
      
      if(any(!is.na(pSet@sensitivity$raw))) {
        if(!all(dimnames(pSet@sensitivity$raw)[[1]] %in% rownames(pSet@sensitivity$info))) {
          warning("For some experiments there is raw sensitivity data but no experimet information in sensitivity info")
        }
      }
      if(!all(rownames(pSet@sensitivity$profiles) %in% rownames(pSet@sensitivity$info))) {
        warning("For some experiments there is sensitivity profiles but no experimet information in sensitivity info")
      }
    }
  }

