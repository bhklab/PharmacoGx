# ==== PharmacoSet Class

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
#' CCLEdrugs  <- treatmentNames(CCLEsmall)
#' CCLEcells <- sampleNames(CCLEsmall)
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
setMethod('subsetTo', signature(object='PharmacoSet'), function(object,
  cells=NULL, drugs=NULL, molecular.data.cells=NULL, keep.controls=TRUE, ...)
{
    .subsetToPharmacoSet(object, cells=cells, drugs=drugs,
        molecular.data.cells=molecular.data.cells, keep.controls=keep.controls)
})


#' @importFrom CoreGx .intersectList
#' @keywords internal
.subsetToPharmacoSet <- function(object, cells=NULL, drugs=NULL,
    molecular.data.cells=NULL, keep.controls=TRUE, ...)
{
  drop=FALSE #TODO:: Is this supposed to be here?

  adArgs = list(...)
  if ('exps' %in% names(adArgs)) {
  	exps <- adArgs[['exps']]
  	if(is(exps, 'data.frame')) {
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

  molecularProfilesSlot(object) <- lapply(molecularProfilesSlot(object), function(SE, cells, drugs, molecular.data.cells){

    molecular.data.type <- ifelse(length(grep('rna', S4Vectors::metadata(SE)$annotation) > 0), 'rna', S4Vectors::metadata(SE)$annotation)
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cells <- molecular.data.cells[[molecular.data.type]]
    }

        column_indices <- NULL

      if (length(cells)==0 && length(drugs)==0) {
          column_indices <- seq_len(ncol(SE)) # This still returns the number of samples in an SE, but without a label
      }
      if(length(cells)==0 && datasetType(object)=='sensitivity') {
        column_indices <- seq_len(ncol(SE))
      }

      cell_line_index <- NULL
      if(length(cells)!=0) {
        if (!all(cells %in% sampleNames(object))) {
              stop('Some of the cell names passed to function did not match to names in the PharmacoSet. Please ensure you are using cell names as returned by the cellNames function')
        }
          cell_line_index <- which(SummarizedExperiment::colData(SE)[["sampleid"]] %in% cells)
        # if (length(na.omit(cell_line_index))==0){
    #       stop('No cell lines matched')
    #     }
      }
      drugs_index <- NULL
      if(datasetType(object)=='perturbation' || datasetType(object)=='both'){
        if(length(drugs) != 0) {
            if (!all(drugs %in% treatmentNames(object))){
                  stop('Some of the drug names passed to function did not match to names in the PharmacoSet. Please ensure you are using drug names as returned by the drugNames function')
            }
          drugs_index <- which(SummarizedExperiment::colData(SE)[["treatmentid"]] %in% drugs)
          # if (length(drugs_index)==0){
    #         stop('No drugs matched')
    #       }
          if(keep.controls) {
            control_indices <- which(SummarizedExperiment::colData(SE)[['xptype']]=='control')
            drugs_index <- c(drugs_index, control_indices)
          }
        }
      }

      if(length(drugs_index) != 0 && length(cell_line_index) != 0) {
        if(length(intersect(drugs_index, cell_line_index)) == 0) {
          stop('This Drug - Cell Line combination was not tested together.')
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

  if ((datasetType(object) == 'sensitivity' | datasetType(object) == 'both') & length(exps) != 0) {
      sensitivityInfo(object) <- sensitivityInfo(object)[exps, , drop=drop]
      rownames(sensitivityInfo(object)) <- names(exps)
      if(length(sensitivityRaw(object)) > 0) {
        sensitivityRaw(object) <- sensitivityRaw(object)[exps, , , drop=drop]
        dimnames(sensitivityRaw(object))[[1]] <- names(exps)
      }
      sensitivityProfiles(object) <- sensitivityProfiles(object)[exps, , drop=drop]
      rownames(sensitivityProfiles(object)) <- names(exps)

      sensNumber(object) <- .summarizeSensitivityNumbers(object)
  }
  else if ((datasetType(object) == 'sensitivity' | datasetType(object) == 'both') & (length(drugs) != 0 | length(cells) != 0)) {

        drugs_index <- which(sensitivityInfo(object)[, "treatmentid"] %in% drugs)
        cell_line_index <- which(sensitivityInfo(object)[,"sampleid"] %in% cells)
        if (length(drugs_index) !=0 & length(cell_line_index) !=0 ) {
          if (length(intersect(drugs_index, cell_line_index)) == 0) {
            stop('This Drug - Cell Line combination was not tested together.')
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
        treatmentResponse(object)[names(treatmentResponse(object))[names(treatmentResponse(object))!='n']] <- lapply(treatmentResponse(object)[names(treatmentResponse(object))[names(treatmentResponse(object))!='n']], function(x,i, drop){
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
		if(datasetType(object) == 'sensitivity' | datasetType(object) == 'both'){
			drugs <- unique(sensitivityInfo(object)[["treatmentid"]])
		}
		if(datasetType(object) == 'perturbation' | datasetType(object) == 'both'){
			drugs <- union(drugs, na.omit(CoreGx::.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(colData(SE)[["treatmentid"]])}))))
		}
	}
	if (length(cells)==0) {
		cells <- union(cells, na.omit(CoreGx::.unionList(lapply(molecularProfilesSlot(object), function(SE){unique(colData(SE)[["sampleid"]])}))))
        if (datasetType(object) =='sensitivity' | datasetType(object) == 'both'){
            cells <- union(cells, sensitivityInfo(object)[["sampleid"]])
        }
	}
	treatmentInfo(object) <- treatmentInfo(object)[drugs, , drop=drop]
	sampleInfo(object) <- sampleInfo(object)[cells, , drop=drop]
	curation(object)$treatment <- curation(object)$treatment[drugs, , drop=drop]
	curation(object)$sample <- curation(object)$sample[cells, , drop=drop]
	curation(object)$tissue <- curation(object)$tissue[cells, , drop=drop]
	if (datasetType(object) == 'sensitivity' | datasetType(object) == 'both'  & length(exps) == 0) {
	  sensNumber(object) <- sensNumber(object)[cells, drugs, drop=drop]
	}
	if (datasetType(object) == 'perturbation' | datasetType(object) == 'both') {
	  pertNumber(object) <- pertNumber(object)[cells, drugs, , drop=drop]
    }
  return(object)
}