#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export
NULL

setClassUnion('list_OR_MAE', c('list', 'MultiAssayExperiment'))

# #' @importClassesFrom CoreGx LongTable TreatmentResponseExperiment
# setClassUnion('list_OR_LongTable', c('list', 'LongTable'))

.local_class="PharmacoSet"

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
#' @slot sample A \code{data.frame} containing the annotations for all the cell
#'   lines profiled in the data set, across all data types
#' @slot treatment A \code{data.frame} containg the annotations for all the drugs
#'   profiled in the data set, across all data types
#' @slot treatmentResponse A \code{list} containing all the data for the
#'   sensitivity experiments, including \code{$info}, a \code{data.frame}
#'   containing the experimental info,\code{$raw} a 3D \code{array} containing
#'   raw data, \code{$profiles}, a \code{data.frame} containing sensitivity
#'   profiles statistics, and \code{$n}, a \code{data.frame} detailing the
#'   number of experiments for each cell-drug pair
#' @slot perturbation A \code{list} containting \code{$n}, a \code{data.frame}
#'   summarizing the available perturbation data,
#' @slot curation A \code{list} containing mappings for \code{$treatment},
#'   \code{cell}, \code{tissue} names  used in the data set to universal
#'   identifiers used between different PharmacoSet objects
#' @slot datasetType A \code{character} string of 'sensitivity',
#'   'perturbation', or both detailing what type of data can be found in the
#'   PharmacoSet, for proper processing of the data
#'
#' @importClassesFrom CoreGx CoreSet
#' @importClassesFrom CoreGx LongTable
#' @importClassesFrom CoreGx TreatmentResponseExperiment
#'
#' @return An object of the PharmacoSet class
.PharmacoSet <- setClass('PharmacoSet',
    contains='CoreSet')


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
#' @inheritParams CoreGx::CoreSet
#'
#' @return An object of class `PharmacoSet`
#
#' @import methods
#' @importFrom utils sessionInfo
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment rowData colData assay assays assayNames Assays
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom CoreGx CoreSet
#'
#' @export
PharmacoSet <-  function(name, molecularProfiles=list(), sample=data.frame(),
        treatment=data.frame(), sensitivityInfo=data.frame(),
        sensitivityRaw=array(dim=c(0,0,0)), sensitivityProfiles=matrix(),
        sensitivityN=matrix(nrow=0, ncol=0), perturbationN=array(NA, dim=c(0,0,0)),
        curationTreatment=data.frame(), curationSample = data.frame(),
        curationTissue = data.frame(), datasetType=c("sensitivity", "perturbation", "both"),
        verify = TRUE, ...) {

    #.Deprecated("PharmacoSet2", )

    cSet <- CoreGx::CoreSet(
        name=name,
        molecularProfiles = molecularProfiles,
        sample=sample,
        treatment=treatment,
        sensitivityInfo=sensitivityInfo,
        sensitivityRaw=sensitivityRaw,
        sensitivityProfiles=sensitivityProfiles,
        sensitivityN=sensitivityN,
        perturbationN=perturbationN,
        curationTreatment=curationTreatment,
        curationSample=curationSample,
        curationTissue=curationTissue,
        datasetType=datasetType,
        verify=verify,
        ...
    )

    pSet  <- .PharmacoSet(
        annotation=cSet@annotation,
        molecularProfiles=cSet@molecularProfiles,
        sample=cSet@sample,
        treatment=cSet@treatment,
        datasetType=cSet@datasetType,
        treatmentResponse=cSet@treatmentResponse,
        perturbation=cSet@perturbation,
        curation=cSet@curation
    )
    if (verify) checkPsetStructure(pSet)
    if (length(sensitivityN) == 0 && datasetType %in% c("sensitivity", "both")) {
        pSet@treatmentResponse$n <- .summarizeSensitivityNumbers(pSet)
    }
    if (!length(perturbationN) &&
            datasetType %in% c("perturbation", "both")) {
        pSet@perturbation$n <- .summarizePerturbationNumbers(pSet)
    }
    return(pSet)
}

#' @eval CoreGx:::.docs_CoreSet2_constructor(class_=.local_class,
#' sx_="Samples in a `PharmacoSet` represent cancer cell-lines.",
#' tx_="Treatments in a `PharmacoSet` represent pharmaceutical compounds.",
#' cx_="This class requires an additional curation item, tissue, which maps
#' from published to standardized tissue idenifiers.",
#' data_=.local_data)
#' @importFrom CoreGx CoreSet2 LongTable TreatmentResponseExperiment
#' @export
PharmacoSet2 <- function(name="emptySet", treatment=data.frame(),
        sample=data.frame(), molecularProfiles=MultiAssayExperiment(),
        treatmentResponse=TreatmentResponseExperiment(),
        perturbation=list(),
        curation=list(sample=data.frame(), treatment=data.frame(),
        tissue=data.frame()), datasetType="sensitivity"
) {
    # -- Leverage existing checks in CoreSet constructor
    cSet <- CoreSet2(name=name, treatment=treatment,
        sample=sample, treatmentResponse=treatmentResponse,
        molecularProfiles=molecularProfiles, curation=curation,
        perturbation=perturbation, datasetType=datasetType)

    ## -- data integrity
    # treatment
    ## TODO

    .PharmacoSet(
        annotation=cSet@annotation,
        sample=cSet@sample,
        treatment=cSet@treatment,
        molecularProfiles=cSet@molecularProfiles,
        treatmentResponse=cSet@treatmentResponse,
        datasetType=cSet@datasetType,
        curation=cSet@curation,
        perturbation=cSet@perturbation
    )
}

# Constructor Helper Functions ----------------------------------------------

#' @keywords internal
#' @importFrom CoreGx idCols . .errorMsg .collapse
.summarizeSensitivityNumbers <- function(object) {
    ## TODO:: Checks don't like assigning to global evnironment. Can we return this?
    assign('object_sumSenNum', object) # Removed envir=.GlobalEnv
    if (datasetType(object) != 'sensitivity' && datasetType(object) != 'both') {
        stop ('Data type must be either sensitivity or both')
    }
    ## consider all drugs
    drugn <- treatmentNames(object)
    ## consider all cell lines
    celln <- sampleNames(object)
    sensitivity.info <- matrix(0, nrow=length(celln), ncol=length(drugn),
        dimnames=list(celln, drugn))
    drugids <- sensitivityInfo(object)[ , "treatmentid"]
    sampleids <- sensitivityInfo(object)[ , "sampleid"]
    sampleids <- sampleids[grep('///', drugids, invert=TRUE)]
    drugids <- drugids[grep('///', drugids, invert=TRUE)]
    tt <- table(sampleids, drugids)
    sensitivity.info[rownames(tt), colnames(tt)] <- tt

    return(sensitivity.info)
}

#' @importFrom CoreGx .summarizeMolecularNumbers
.summarizeMolecularNumbers <- function(object) {
    CoreGx::.summarizeMolecularNumbers(object)
}

#' @importFrom CoreGx treatmentNames sampleNames
.summarizePerturbationNumbers <- function(object) {

    if (datasetType(object) != 'perturbation' && datasetType(object) != 'both') {
        stop('Data type must be either perturbation or both')
    }

    ## consider all drugs
    drugn <- treatmentNames(object)

    ## consider all cell lines
    celln <- sampleNames(object)

    mprof <- molecularProfilesSlot(object)
    perturbation.info <- array(0, dim=c(length(celln), length(drugn),
        length(mprof)),
        dimnames=list(celln, drugn, names(mprof))
    )
    for (i in seq_len(length(mprof))) {
        if (nrow(colData(mprof[[i]])) > 0 &&
                all(c("sampleid", "treatmentid") %in%
                    colnames(mprof[[i]]))) {
            tt <- table(
                colData(mprof[[i]])[, "sampleid"],
                colData(mprof[[i]])[, "treatmentid"]
            )
        perturbation.info[rownames(tt), colnames(tt), names(mprof)[i]] <- tt
        }
    }

    return(perturbation.info)
}

### -------------------------------------------------------------------------
### Class Validity ----------------------------------------------------------
### -------------------------------------------------------------------------

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
  function(object, plotDist=FALSE, result.dir='.') {

    # Make directory to store results if it doesn't exist
    if(!file.exists(result.dir) & plotDist) { dir.create(result.dir, showWarnings=FALSE, recursive=TRUE) }

    #####
    # Checking molecularProfiles
    #####
    # Can this be parallelized or does it mess with the order of printing warnings?
    mprof <- molecularProfilesSlot(object)
    for( i in seq_along(mprof)) {
      profile <- mprof[[i]]
      nn <- names(mprof)[i]

      # Testing plot rendering for rna and rnaseq
      if((S4Vectors::metadata(profile)$annotation == 'rna' || S4Vectors::metadata(profile)$annotation == 'rnaseq') && plotDist)
      {
        pdf(file=file.path(result.dir, sprintf('%s.pdf', nn)))
        hist(assays(profile)[[1]], breaks = 100)
        dev.off()
      }

      ## Test if sample and feature annotations dimensions match the assay
      warning(ifelse(nrow(rowData(profile)) != nrow(assays(profile)[[1]]),
                     sprintf('%s: number of features in fData is different from
                             SummarizedExperiment slots', nn),
                     sprintf('%s: rowData dimension is OK', nn)
                     )
              )
      warning(ifelse(nrow(colData(profile)) != ncol(assays(profile)[[1]]),
                     sprintf('%s: number of cell lines in pData is different
                             from expression slots', nn),
                     sprintf('%s: colData dimension is OK', nn)
                     )
              )


      # Checking sample metadata for required columns
      warning(ifelse("sampleid" %in% colnames(colData(profile)), '',
                     sprintf('%s: sampleid does not exist in colData (samples)
                             columns', nn)))
      warning(ifelse('batchid' %in% colnames(colData(profile)), '',
                     sprintf('%s: batchid does not exist in colData (samples)
                             columns', nn)))

      # Checking mDataType of the SummarizedExperiment for required columns
      if(S4Vectors::metadata(profile)$annotation == 'rna' |
         S4Vectors::metadata(profile)$annotation == 'rnaseq')
      {
        warning(ifelse('BEST' %in% colnames(rowData(profile)), 'BEST is OK',
                       sprintf('%s: BEST does not exist in rowData (features)
                               columns', nn)))
        warning(ifelse('Symbol' %in% colnames(rowData(profile)), 'Symbol is OK',
                       sprintf('%s: Symbol does not exist in rowData (features)
                               columns', nn)))
      }

      # Check that all sampleids from the object are included in molecularProfiles
      if("sampleid" %in% colnames(colData(profile))) {
        if (!all(colData(profile)[,"sampleid"] %in% sampleNames(object))) {
          warning(sprintf('%s: not all the cell lines in this profile are in
                          cell lines slot', nn))
        }
      }else {
        warning(sprintf('%s: sampleid does not exist in colData (samples)', nn))
      }
    }

}


### -------------------------------------------------------------------------
### Method Definitions ------------------------------------------------------
### -------------------------------------------------------------------------

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
setMethod('show', signature=signature(object='PharmacoSet'), function(object) {
    callNextMethod(object)
})

#' Get the dimensions of a PharmacoSet
#'
#' @param x PharmacoSet
#' @return A named vector with the number of Cells and Drugs in the PharmacoSet
#' @export
setMethod('dim', signature=signature(x='PharmacoSet'), function(x){
    return(c(Cells=length(sampleNames(x)), Drugs=length(treatmentNames(x))))
})


### TODO:: Add updating of sensitivity Number tables
#' @importFrom CoreGx updateSampleId
#' @aliases updateCellId
updateSampleId <- updateCellId <- function(object, new.ids = vector('character')){
    CoreGx::updateSampleId(object, new.ids)
}
