########################
## Petr Smirnov
## All rights Reserved
## September 1, 2014
########################

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
#' 
#' @param pSet A \code{PharmacoSet} object
#' @param object A \code{PharmacoSet} object
#' @param value A replacement value
#' 
#' @slot annotation A \code{list} of annotation data about the PharmacoSet,
#'    including the \code{$name} and the session information for how the object
#'    was creating, detailing the exact versions of R and all the packages used
#' @slot molecularData A \code{list} containing 4 \code{Biobase::ExpressionSet} 
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
#'   \code{$phenotype}, a \code{data.frame} containing sensitivity phenotype 
#'   statistics, and \code{$n}, a \code{data.frame} detailing the number of 
#'   experiments for each cell-drug pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame} 
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for \code{$drug}, 
#'   \code{cell}, \code{tissue} names  used in the data set to universal 
#'   identifiers used between different PharmacoSet objects
#' @slot datasetType A \code{character} string of 'sensitivity', 
#'   'preturbation', or both detailing what type of data can be found in the 
#'   PharmacoSet, for proper processing of the data
.PharmacoSet <- setClass("PharmacoSet", slots = list(
													 annotation = "list",
													 molecularData = "list",
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
#test <- new("PharmacoSet",cell=data.cmap$cell, drug=data.cmap$drug)
#exprs(test)
# setGeneric("PharmacoSet",
# 	function(annotation,
# 		molecularData,
# 		cell=matrix(),
# 		drug=list(),
# 		sensitivity=list(),
# 		datasetType=c("sensitivity", "perturbation"),
# 		... ) standardGeneric("PharmacoSet"), signature="ExpressionSet")

#' PharmacoSet constructor
#' 
#' A constructor that simplifies the process of creating PharmacoSets, as well 
#' as creates empty objects for data not provided to the constructor. Only
#' objects returned by this constructor are expected to work with the PharmacoSet
#' methods.
#' 
#' @param name A \code{character} string detailing the name of the dataset
#' @param rna,dna,snp,cnv A \code{Biobase::ExpressionSet} type object containing
#'   data for RNA, DNA SNP and Copy Number Variation respectively, with
#'   associated \code{fData} and \code{pData} containing the row and column
#'   metadata
#' @param cell A \code{data.frame} containg the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @param drug A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @param sensitivityInfo A \code{data.frame} containing the information for the
#'   sensitivity experiments
#' @param sensitivityRaw A 3 Dimensional \code{array} contaning the raw drug
#'   dose response data for the sensitivity experiments
#' @param sensitivityPhenotype \code{data.frame} containing drug sensitivity phenotype 
#'   statistics such as IC50 and AUC
#' @param sensitivityN,perturbationN A \code{data.frame} summarizing the
#'   available sensitivity/perturbation data
#' @param curationDrug,curationCell,curationTissue A \code{data.frame} mapping
#'   the names for drugs, cells and tissues used in the data set to universal
#'   identifiers used between different PharmacoSet objects
#' @param datasetType A \code{character} string of 'sensitivity',
#'   'preturbation', or both detailing what type of data can be found in the
#'   PharmacoSet, for proper processing of the data
#'
#' @return An object of class PharmacoSet
#' @export
#' @import Biobase
#' @import methods
PharmacoSet<-  function(name, 
                          rna=ExpressionSet(), 
                          dna=ExpressionSet(), 
                          snp=ExpressionSet(), 
                          cnv=ExpressionSet(), 
                          cell=data.frame(), 
                          drug=data.frame(), 
                          sensitivityInfo=data.frame(),
                          sensitivityRaw=array(dim=c(0,0,0)), 
                          sensitivityPhenotype=matrix(), 
                          sensitivityN=matrix(nrow=0, ncol=0), 
                          perturbationN=array(NA, dim=c(0,0,0)), 
                          curationDrug=data.frame(), 
                          curationCell = data.frame(), 
                          curationTissue = data.frame(), 
                          datasetType=c("sensitivity", "perturbation", "both"))
{
	datasetType <- match.arg(datasetType)
	
	annotation <- list()
	annotation$name <- as.character(name)
	annotation$dateCreated <- date()
	annotation$sessionInfo <- sessionInfo()
	annotation$call <- match.call()
	
	molecularData <- list("dna"=dna, "rna"=rna, "snp"=snp, "cnv"=cnv)
	for (i in 1:length(molecularData)){
		if (class(molecularData[[i]]) != "ExpressionSet"){
			stop(sprintf("Please provide the %s data as an ExpressionSet", names(molecularData[i])))
		}
	}
	#if (class(cell)!="data.frame"){
	#	stop("Please provide the cell line annotations as a data frame.")
	#}
	#if (class(drug)!="data.frame"){
	#	stop("Please provide the drug annotations as a data frame.")
	#}
	
	sensitivity <- list()
	
	if (!all(rownames(sensitivityInfo) == rownames(sensitivityPhenotype) & rownames(sensitivityInfo) == dimnames(sensitivityRaw)[[1]])){
		stop("Please ensure all the row names match between the sensitivity data.")
	}
	
	sensitivity$info <- sensitivityInfo
	sensitivity$raw <- sensitivityRaw
	sensitivity$phenotype <- sensitivityPhenotype
	sensitivity$n <- sensitivityN
	
	curation <- list()
	curation$drug <- curationDrug
	curation$cell <- curationCell
	curation$tissue <- curationTissue
	### TODO:: Make sure to fix the curation to check for matching row names to the drug and cell line matrices!!!!!!
	
	
	perturbation <- list()
	perturbation$n <- perturbationN
	if (datasetType == "perturbation" || datasetType == "both") {
		perturbation$info <- "The metadata for the perturbation experiments is available for each molecular type by calling the appropriate info function. \n For example, for RNA transcriptome perturbations, the metadata can be accessed using rnaInfo(pSet)."
	} else {
		perturbation$info <- "Not a perturbation dataset."
	}
	
	pSet  <- .PharmacoSet(annotation=annotation, molecularData=molecularData, cell=as.data.frame(cell), drug=as.data.frame(drug), datasetType=datasetType, sensitivity=sensitivity, perturbation=perturbation, curation=curation)
	return(pSet)
}
	
#' cellInfo Generic
#' 
#' Generic for cellInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve cell info from
#' @return a \code{data.frame} with the cell annotations
setGeneric("cellInfo", function(pSet) standardGeneric("cellInfo"))

#' cellInfo<- Generic
#' 
#' Generic for cellInfo replace method
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
#' @param pSet The \code{PharmacoSet} to retrieve drug info from
#' @return a \code{data.frame} with the drug annotations
setGeneric("drugInfo", function(pSet) standardGeneric("drugInfo"))

#' drugInfo<- Generic
#' 
#' Generic for drugInfo replace method
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

#' rnaInfo Generic
#' 
#' Generic for rnaInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve rna annotations from
#' @return a \code{data.frame} with the rna annotations
setGeneric("rnaInfo", function(pSet) standardGeneric("rnaInfo"))

#' @describeIn PharmacoSet Return the RNA experiment info from the PharmacoSet 
#' @export
setMethod(rnaInfo, "PharmacoSet", function(pSet){
	
	return(pData(pSet@molecularData$rna))
	
})

#' rnaInfo<- Generic
#' 
#' Generic for rnaInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace rna annotations in
#' @param value A \code{data.frame} with the new rna annotations
#' @return Updated \code{PharmacoSet}
setGeneric("rnaInfo<-", function(object, value) standardGeneric("rnaInfo<-"))
#' @describeIn PharmacoSet Update the RNA experiment info from the PharmacoSet 
#' @export
setReplaceMethod("rnaInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

	pData(object@molecularData$rna) <- value
	object
})

#' rnaData Generic
#' 
#' Generic for rnaData method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve rna data from
#' @return a \code{data.frame} with the rna experiment data
setGeneric("rnaData", function(pSet) standardGeneric("rnaData"))
#' @describeIn PharmacoSet Return the RNA data from the PharmacoSet 
#' @export
setMethod(rnaData, "PharmacoSet", function(pSet){
	
	return(exprs(pSet@molecularData$rna))
	
})

#' rnaData<- Generic
#' 
#' Generic for rnaData replace method
#' 
#' @param object The \code{PharmacoSet} to replace rna data in
#' @param value A \code{matrix} with the new rna data
setGeneric("rnaData<-", function(object, value) standardGeneric("rnaData<-"))
#' @describeIn PharmacoSet Update the RNA data from the PharmacoSet 
#' @export
setReplaceMethod("rnaData", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

	exprs(object@molecularData$rna) <- value
	object
})

#' dnaInfo Generic
#' 
#' Generic for dnaInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve dna annotations from
#' @return a \code{data.frame} with the dna annotations
setGeneric("dnaInfo", function(pSet) standardGeneric("dnaInfo"))
#' @describeIn PharmacoSet Return the DNA experiment info from the PharmacoSet 
#' @export
setMethod(dnaInfo, "PharmacoSet", function(pSet){
	
	return(pData(pSet@molecularData$dna))
	
})

#' dnaInfo<- Generic
#' 
#' Generic for dnaInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace dna info in
#' @param value A \code{data.frame} with the new dna annotations
#' @return Updated \code{PharmacoSet}
setGeneric("dnaInfo<-", function(object, value) standardGeneric("dnaInfo<-"))

#' @describeIn PharmacoSet UPdate the DNA experiment info from the PharmacoSet 
#' @export
setReplaceMethod("dnaInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

	pData(object@molecularData$dna) <- value
	object
})

#' dnaData Generic
#' 
#' Generic for dnaData method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve dna data from
#' @return a \code{data.frame} with the dna experiment data
setGeneric("dnaData", function(pSet) standardGeneric("dnaData"))
#' @describeIn PharmacoSet Return the DNA data from the PharmacoSet 
#' @export
setMethod(dnaData, "PharmacoSet", function(pSet){
	
	return(exprs(pSet@molecularData$dna))
	
})

#' dnaData<- Generic
#' 
#' Generic for dnaData replace method
#' 
#' @param object The \code{PharmacoSet} to replace dna data in
#' @param value A \code{matrix} with the new dna data
setGeneric("dnaData<-", function(object, value) standardGeneric("dnaData<-"))

#' @describeIn PharmacoSet Update the DNA data from the PharmacoSet 
#' @export
setReplaceMethod("dnaData", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

	exprs(object@molecularData$dna) <- value
	object
})

#' snpInfo Generic
#' 
#' Generic for snpInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve snp annotations from
#' @return a \code{data.frame} with the snp annotations
setGeneric("snpInfo", function(pSet) standardGeneric("snpInfo"))

#' @describeIn PharmacoSet Return the SNP experiment info from the PharmacoSet 
#' @export
setMethod(snpInfo, "PharmacoSet", function(pSet){
	
	return(pData(pSet@molecularData$snp))
	
})

#' snpInfo<- Generic
#' 
#' Generic for snpInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace snp annotations in
#' @param value A \code{data.frame} with the new snp annotations
#' @return Updated \code{PharmacoSet}
setGeneric("snpInfo<-", function(object, value) standardGeneric("snpInfo<-"))

#' @describeIn PharmacoSet Update the SNP experiment info from the PharmacoSet 
#' @export
setReplaceMethod("snpInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

	pData(object@molecularData$snp) <- value
	object
})

#' snpData Generic
#' 
#' Generic for snpData method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve snp data from
#' @return a \code{data.frame} with the snp experiment data
setGeneric("snpData", function(pSet) standardGeneric("snpData"))

#' @describeIn PharmacoSet Return the SNP data from the PharmacoSet 
#' @export
setMethod(snpData, "PharmacoSet", function(pSet){
	
	return(exprs(pSet@molecularData$snp))
	
})

#' snpData<- Generic
#' 
#' Generic for snpData replace method
#' 
#' @param object The \code{PharmacoSet} to replace snp data in
#' @param value A \code{matrix} with the new snp data
setGeneric("snpData<-", function(object, value) standardGeneric("snpData<-"))

#' @describeIn PharmacoSet Update the SNP data from the PharmacoSet 
#' @export
setReplaceMethod("snpData", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

	exprs(object@molecularData$snp) <- value
	object
})

#' cnvInfo Generic
#' 
#' Generic for cnvInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve cnv annotations from
#' @return a \code{data.frame} with the cnv annotations
setGeneric("cnvInfo", function(pSet) standardGeneric("cnvInfo"))
#' @describeIn PharmacoSet Return the CNV experiment info from the PharmacoSet 
#' @export
setMethod(cnvInfo, "PharmacoSet", function(pSet){
	
	return(pData(pSet@molecularData$cnv))
	
})

#' cnvInfo<- Generic
#' 
#' Generic for cnvInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace cnv annotations in
#' @param value A \code{data.frame} with the new cnv annotations
#' @return Updated \code{PharmacoSet}
setGeneric("cnvInfo<-", function(object, value) standardGeneric("cnvInfo<-"))

#' @describeIn PharmacoSet Update the CNV experiment info from the PharmacoSet 
#' @export
setReplaceMethod("cnvInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

	pData(object@molecularData$cnv) <- value
	object
})

#' cnvData Generic
#' 
#' Generic for cnvData method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve cnv data from
#' @return a \code{data.frame} with the cnv experiment data
setGeneric("cnvData", function(pSet) standardGeneric("cnvData"))
#' @describeIn PharmacoSet Return the CNV data from the PharmacoSet 
#' @export
setMethod(cnvData, "PharmacoSet", function(pSet){
	
	return(exprs(pSet@molecularData$cnv))
	
})

#' cnvData<- Generic
#' 
#' Generic for cnvData replace method
#' 
#' @param object The \code{PharmacoSet} to replace cnv data in
#' @param value A \code{matrix} with the new cnv data
setGeneric("cnvData<-", function(object, value) standardGeneric("cnvData<-"))

#' @describeIn PharmacoSet Update the CNV data from the PharmacoSet 
#' @export
setReplaceMethod("cnvData", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

	exprs(object@molecularData$cnv) <- value
	object
})

#' geneInfo Generic
#' 
#' Generic for geneInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve gene annotations from
#' @return a \code{data.frame} with the gene annotations
setGeneric("geneInfo", function(pSet) standardGeneric("geneInfo"))

#' @describeIn PharmacoSet Return the gene info for the molecular data 
#' @export
setMethod(geneInfo, "PharmacoSet", function(pSet){
	#### Do this for all datasets????????????
	return(fData(pSet@molecularData$rna))
	
})

#' geneInfo<- Generic
#' 
#' Generic for geneInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace gene annotations in
#' @param value A \code{data.frame} with the new gene annotations
#' @return Updated \code{PharmacoSet}
setGeneric("geneInfo<-", function(object, value) standardGeneric("geneInfo<-"))

#' @describeIn PharmacoSet Replace the gene info for the molecular data
#' @export
setReplaceMethod("geneInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){
	
	fData(object@molecularData$rna) <- value
	
	object
})

#' sensitivityInfo Generic
#' 
#' Generic for rnaInfo method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve sensitivity experiment annotations from
#' @return a \code{data.frame} with the sensitivity annotations
setGeneric("sensitivityInfo", function(pSet) standardGeneric("sensitivityInfo"))
#' @describeIn PharmacoSet Return the sensitivity experiment info
#' @export
setMethod(sensitivityInfo, "PharmacoSet", function(pSet){
	
	return(pSet@sensitivity$info)
	
})

#' sensitivityInfo<- Generic
#' 
#' Generic for sensitivityInfo replace method
#' 
#' @param object The \code{PharmacoSet} to replace sensitivity experiment annotations in
#' @param value A \code{data.frame} with the new sensitivity annotations
#' @return Updated \code{PharmacoSet}
setGeneric("sensitivityInfo<-", function(object, value) standardGeneric("sensitivityInfo<-"))

#' @describeIn PharmacoSet Update the sensitivity experiment info
#' @export
setReplaceMethod("sensitivityInfo", signature = signature(object="PharmacoSet",value="data.frame"), function(object, value){

	object@sensitivity$info <- value
	object
})

#' sensitivityData Generic
#' 
#' Generic for sensitivityData method 
#' 
#' @param pSet The \code{PharmacoSet} to retrieve sensitivity data from
#' @return a \code{data.frame} with the sensitivity experiment data
setGeneric("sensitivityData", function(pSet) standardGeneric("sensitivityData"))
#' @describeIn PharmacoSet Return the phenotypic data for the drug dose sensitivity
#' @export
setMethod(sensitivityData, "PharmacoSet", function(pSet){
	
	return(pSet@sensitivity$phenotype)
	
})

#' sensitivityData<- Generic
#' 
#' Generic for sensitivityData replace method
#' 
#' @param object The \code{PharmacoSet} to replace sensitivity data in
#' @param value A \code{matrix} with the new sensitivity data
setGeneric("sensitivityData<-", function(object, value) standardGeneric("sensitivityData<-"))

#' @describeIn PharmacoSet Update the phenotypic data for the drug dose
#'   sensitivity
#' @export
setReplaceMethod("sensitivityData", signature = signature(object="PharmacoSet",value="matrix"), function(object, value){

	object@sensitivity$phenotype <- value
	object
})

#' sensitivityMeasures Generic
#' 
#' A generic for the sensitivity measures function
#'
#' @param pSet A \code{PharmacoSet}
#' @return A \code{character} vector of available sensitivity measures 
setGeneric("sensitivityMeasures", function(pSet) standardGeneric("sensitivityMeasures"))
#' @describeIn PharmacoSet Returns the available sensitivity phenotype
#'   summaries, for example, whether there are IC50 values available
#' @export
setMethod(sensitivityMeasures, "PharmacoSet", function(pSet){
	
	return(colnames(sensitivityData(pSet)))
	
})

#' drugNames Generic
#' 
#' A generic for the drugNames method
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
  rownames(pSet@drug)

})



#' drugNames<- Generic
#' 
#' A generic for the drugNames replacement method
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
#' @param pSet The \code{PharmacoSet} to return cell names from
#' @return A vector of the cell names used in the PharmacoSet
setGeneric("cellNames", function(pSet) standardGeneric("cellNames"))
#' @describeIn PharmacoSet Return the cell names used in the dataset
#' @export
setMethod(cellNames, "PharmacoSet", function(pSet){
  
  rownames(pSet@cell)
  
})

#' cellNames<- Generic
#' 
#' A generic for the cellNames replacement method
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

#' geneNames Generic
#' 
#' A generic for the geneNames method
#' 
#' @param pSet The \code{PharmacoSet} to return gene names from
#' @return A vector of the gene names used in the PharmacoSet
setGeneric("geneNames", function(pSet) standardGeneric("geneNames"))

#' @describeIn PharmacoSet Return the gene names used in the dataset
#' @export
setMethod(geneNames, "PharmacoSet", function(pSet){
  
  rownames(geneInfo(pSet))
  
})

#' dateCreated Generic
#' 
#' A generic for the dateCreated method
#' 
#' @param pSet A \code{PharmacoSet} 
#' @return The date the PharmacoSet was created
setGeneric("dateCreated", function(pSet) standardGeneric("dateCreated"))
#' @describeIn PharmacoSet Return the date the PharmacoSet was created
#' @export
setMethod(dateCreated, "PharmacoSet", function(pSet) {
  PharmacoSet@annotation$dateCreated
})

#' pSetName Generic
#' 
#' A generic for the pSetName method
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
#' @export
setMethod("show", signature=signature(object="PharmacoSet"), 
	function(object) {
		cat("Name: ", object@annotation$name, "\n")
		cat("Date Created: ", object@annotation$dateCreated, "\n")
    cat("Number of cell lines: ", nrow(object@cell), "\n")
    cat("Number of drug compounds: ", nrow(object@drug), "\n")
		cat("DNA: \n")
		cat("\tDim: ", dim(object@molecularData$dna), "\n")
		cat("RNA: \n")
		cat("\tDim: ", dim(object@molecularData$rna), "\n")
		cat("SNP: \n")
		cat("\tDim: ", dim(object@molecularData$snp), "\n")
		cat("CNV: \n")
		cat("\tDim: ", dim(object@molecularData$cnv), "\n")
		cat("Drug pertubation: \n")
		cat("\tPlease look at pertInfo(pSet) to determine number of experiments for each drug-cell combination.\n")
		cat("Drug sensitivity: \n")
		cat("\tNumber of Experiments: ",nrow(object@sensitivity$raw),"\n")
		cat("\tPlease look at sensInfo(pSet) to determine number of experiments for each drug-cell combination.\n")
	})



#' `[` 
#'  @param x numeric
#'  @param i numeric
#'  @param j numeric
#'  @param ... further arguments
#'  @param drop A boolean flag of whether to drop single dimensions or not
#'  @export
setMethod(`[`, "PharmacoSet", function(x, i, j, ..., drop = FALSE){
	stop("Please use the subsetTo function and the getter functions to access the object data. Other methods of accessing the object slots are discouraged and may lead to unexpected behaviour.")
})

## FIXED? TODO:: Subset function breaks if it doesnt find cell line in sensitivity info
#' A function to subset a PharmacoSet to data containing only specified drugs, cells and genes
#' 
#' This is the prefered method of subsetting a PharmacoSet. This function allows
#' abstraction of the data to the level of biologically relevant objects: drugs,
#' genes and cells. The function will automatically go through all of the
#' combined data in the PharmacoSet and ensure only the requested drugs, genes
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
#' @param genes A list or vector of gene names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all genes will be left in
#'   the dataset.
#' @param keep.controls If the dataset has perturbation type experiments, should
#'   the controls be kept in the dataset? Defaults to true.
#' @return A PharmacoSet with only the selected drugs, cells and genes
#' @export
subsetTo <- function(pSet, cells=NULL, drugs=NULL, genes=NULL, keep.controls=TRUE){
  drop=FALSE
	### TODO:: implement strict subsetting at this level!!!!
  
	### the function missing does not work as expected in the context below, because the arguments are passed to the anonymous
	### function in lapply, so it does not recognize them as missing
  
  pSet@molecularData <- lapply(pSet@molecularData, function(eset, cells, drugs, genes){
	  
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
	  	cell_line_index <- which(pData(eset)[["cellid"]] %in% cells)
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
	      drugs_index <- which(pData(eset)[["drugid"]] %in% drugs)
	      # if (length(drugs_index)==0){
	#         stop("No drugs matched")
	#       }
	      if(keep.controls) {
	        control_indices <- which(pData(eset)[["xptype"]]=="control")
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
  
	  row_indices <- 0:nrow(exprs(eset))
	  if(!length(genes)==0) {
	    row_indices <- which(rownames(fData(eset))%in%genes)
	    # if (length(row_indices)==0){
	#       stop("No genes matched")
	#     }
	  }
  
  
	  eset <- eset[row_indices,column_indices]

	  return(eset)

  }, cells=cells, drugs=drugs, genes=genes)  
  
  if((pSet@datasetType == "sensitivity" | pSet@datasetType == "both")&(length(drugs)!=0|length(cells)!=0)){
	  drugs_index <- which(sensitivityInfo(pSet)[,"drugid"]%in%drugs)
	  cell_line_index <- which(sensitivityInfo(pSet)[,"cellid"]%in%cells)
	  if(length(drugs_index)!=0&length(cell_line_index)!=0){
	    if(length(intersect(drugs_index, cell_line_index))==0){
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
  
	if (length(drugs)==0){
		if(pSet@datasetType == "sensitivity" | pSet@datasetType == "both"){
			drugs <- unique(sensitivityInfo(pSet)[["drugid"]])
		}
		if(pSet@datasetType == "perturbation" | pSet@datasetType == "both"){
			drug <- unionList(drugs, unique(rnaInfo(pSet)[["drugid"]]),unique(dnaInfo(pSet)[["drugid"]]),unique(snpInfo(pSet)[["drugid"]]),unique(cnvInfo(pSet)[["drugid"]]))
		}
	}
	if (length(cells)==0){
		cells <- unionList(unique(sensitivityInfo(pSet)[["cellid"]]), unique(rnaInfo(pSet)[["cellid"]]),unique(dnaInfo(pSet)[["cellid"]]),unique(snpInfo(pSet)[["cellid"]]),unique(cnvInfo(pSet)[["cellid"]]))
	}
	drugInfo(pSet) <- drugInfo(pSet)[drugs,,drop=drop]
	cellInfo(pSet) <- cellInfo(pSet)[cells,,drop=drop]
	pSet@curation$drug <- pSet@curation$drug[drugs,,drop=drop]
	pSet@curation$cell <- pSet@curation$cell[cells,,drop=drop]
	pSet@curation$tissue <- pSet@curation$tissue[cells,,drop=drop]
	if (pSet@datasetType == "sensitivity" | pSet@datasetType == "both") {
	    pSet@sensitivity$n <- pSet@sensitivity$n[cells,drugs,drop=FALSE]
	}
	if(pSet@datasetType == "perturbation" | pSet@datasetType == "both"){
	    pSet@perturbation$n <- pSet@perturbation$n[cells,drugs,,drop=FALSE]
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
  
  
  pSet@molecularData <- lapply(pSet@molecularData, function(eset){
  		
	  myx <- match(pData(eset)[["cellid"]],rownames(cellInfo(pSet)))
	  pData(eset)[["cellid"]]  <- new.ids[myx]
	  return(eset)
  	  })
  myx <- match(rownames(pSet@curation$cell),rownames(cellInfo(pSet)))
  rownames(pSet@curation$cell) <- new.ids[myx]
  rownames(pSet@curation$tissue) <- new.ids[myx]
  if (dim(pertNumber(pSet))[[1]]>0){
    myx <- match(dimnames(pertNumber(pSet))[[1]], rownames(cellInfo(pSet)))
    dimnames(pertNumber(pSet))[[1]] <- new.ids[myx]
  }
  if (nrow(sensNumber(pSet))>0){
    myx <- match(rownames(sensNumber(pSet)), rownames(cellInfo(pSet)))
    rownames(sensNumber(pSet)) <- new.ids[myx]
  }
  rownames(cellInfo(pSet)) <- new.ids
  return(pSet)

}
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
	  pSet@molecularData <- lapply(pSet@molecularData, function(eset){
  		
		  myx <- match(pData(eset)[["drugid"]],rownames(cellInfo(pSet)))
		  pData(eset)[["drugid"]]  <- new.ids[myx]
		  return(eset)
		  })
  }
  myx <- match(rownames(pSet@curation$drug),rownames(drugInfo(pSet)))
  rownames(pSet@curation$drug) <- new.ids[myx]
  if (dim(pertNumber(pSet))[[2]]>0){
    myx <- match(dimnames(pertNumber(pSet))[[2]], rownames(drugInfo(pSet)))
    dimnames(pertNumber(pSet))[[2]] <- new.ids[myx]
  }
  if (ncol(sensNumber(pSet))>0){
    myx <- match(colnames(sensNumber(pSet)), rownames(drugInfo(pSet)))
    colnames(sensNumber(pSet)) <- new.ids[myx]
  }
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
    
  tt <- table(pSet@sensitivity$info[ , "cellid"], pSet@sensitivity$info[ , "drugid"])
  sensitivity.info[rownames(tt), colnames(tt)] <- tt
  
	return(sensitivity.info)
}

.summarizePertubationNumbers <- function(pSet) {

  if (pSet@datasetType != "perturbation" && pSet@datasetType != "both") {
    stop ("Data type must be either perturbation or both")
  }
  
  ## unique drug identifiers
  # drugn <- sort(unique(unlist(lapply(pSet@molecularData, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "drugid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "drugid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all drugs
  drugn <- rownames(pSet@drug)
  
  ## unique cell line identifiers
  # celln <- sort(unique(unlist(lapply(pSet@molecularData, function (x) {
  #   res <- NULL
  #   if (nrow(pData(x)) > 0 & "cellid" %in% colnames(pData(x))) {
  #     res <- pData(x)[ , "cellid"]
  #   }
  #   return (res)
  # }))))
  
  ## consider all cell lines
  celln <- rownames(pSet@cell)
  
  perturbation.info <- array(0, dim=c(length(celln), length(drugn), length(pSet@molecularData)), dimnames=list(celln, drugn, names((pSet@molecularData))))

	for (i in 1:length(pSet@molecularData)) {
	  if (nrow(pData(pSet@molecularData[[i]])) > 0 && all(is.element(c("cellid", "drugid"), colnames(pData(pSet@molecularData[[i]]))))) {
      tt <- table(pData(pSet@molecularData[[i]])[ , "cellid"], pData(pSet@molecularData[[i]])[ , "drugid"])
	    perturbation.info[rownames(tt), colnames(tt), names(pSet@molecularData)[i]] <- tt
	  }
	}
  
	return(perturbation.info)
}
