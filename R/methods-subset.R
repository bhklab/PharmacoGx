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
#' CCLEdrugs  <- drugNames(CCLEsmall)
#' CCLEcells <- cellNames(CCLEsmall)
#' pSet <- subset(CCLEsmall, drugs = CCLEdrugs[1], cells = CCLEcells[1])
#' pSet
#'
#' @param x A \code{PharmacoSet} to be subsetted
#' @param i A list or vector of cell names as used in the dataset to which
#'   the object will be subsetted. If left blank, then all cells will be left in
#'   the dataset.
#' @param j A list or vector of drug names as used in the dataset to which
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
#' @importMethodsFrom BiocGenerics subset
#' @export
setMethod('subset', signature(x='PharmacoSet'),
          function(x, ...) eval(.subsetPharmacoSet(x, ...)))

#function(object, cells, drugs, molecular.data.cells,
#                     keep.controls=TRUE, ...){
#    eval(substitute(.subsetPharmacoSet(object, cells, drugs, molecular.data.cells, keep.controls, ...)))
#})


#' @importFrom CoreGx .intersectList
#' @keywords internal
.subsetPharmacoSet <- function(object, cells=NULL, drugs=NULL, molecular.data.cells=NULL,
                     keep.controls=TRUE, ...) {
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

  object@molecularProfiles <- lapply(object@molecularProfiles, function(SE, cells, drugs, molecular.data.cells){

    molecular.data.type <- ifelse(length(grep('rna', S4Vectors::metadata(SE)$annotation) > 0), 'rna', S4Vectors::metadata(SE)$annotation)
    if (length(grep(molecular.data.type, names(molecular.data.cells))) > 0) {
      cells <- molecular.data.cells[[molecular.data.type]]
    }

        column_indices <- NULL

      if (length(cells)==0 && length(drugs)==0) {
          column_indices <- seq_len(ncol(SE)) # This still returns the number of samples in an SE, but without a label
      }
      if(length(cells)==0 && object@datasetType=='sensitivity') {
        column_indices <- seq_len(ncol(SE))
      }

      cell_line_index <- NULL
      if(length(cells)!=0) {
        if (!all(cells %in% cellNames(object))) {
              stop('Some of the cell names passed to function did not match to names in the PharmacoSet. Please ensure you are using cell names as returned by the cellNames function')
        }
          cell_line_index <- which(SummarizedExperiment::colData(SE)[['cellid']] %in% cells)
        # if (length(na.omit(cell_line_index))==0){
    #       stop('No cell lines matched')
    #     }
      }
      drugs_index <- NULL
      if(object@datasetType=='perturbation' || object@datasetType=='both'){
        if(length(drugs) != 0) {
            if (!all(drugs %in% drugNames(object))){
                  stop('Some of the drug names passed to function did not match to names in the PharmacoSet. Please ensure you are using drug names as returned by the drugNames function')
            }
          drugs_index <- which(SummarizedExperiment::colData(SE)[['drugid']] %in% drugs)
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

  if ((object@datasetType == 'sensitivity' | object@datasetType == 'both') & length(exps) != 0) {
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
  else if ((object@datasetType == 'sensitivity' | object@datasetType == 'both') & (length(drugs) != 0 | length(cells) != 0)) {

        drugs_index <- which (sensitivityInfo(object)[, 'drugid'] %in% drugs)
        cell_line_index <- which (sensitivityInfo(object)[,'cellid'] %in% cells)
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
        object@sensitivity[names(object@sensitivity)[names(object@sensitivity)!='n']] <- lapply(object@sensitivity[names(object@sensitivity)[names(object@sensitivity)!='n']], function(x,i, drop){
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
		if(object@datasetType == 'sensitivity' | object@datasetType == 'both'){
			drugs <- unique(sensitivityInfo(object)[['drugid']])
		}
		if(object@datasetType == 'perturbation' | object@datasetType == 'both'){
			drugs <- union(drugs, na.omit(CoreGx::.unionList(lapply(object@molecularProfiles, function(SE){unique(colData(SE)[['drugid']])}))))
		}
	}
	if (length(cells)==0) {
		cells <- union(cells, na.omit(CoreGx::.unionList(lapply(object@molecularProfiles, function(SE){unique(colData(SE)[['cellid']])}))))
        if (object@datasetType =='sensitivity' | object@datasetType == 'both'){
            cells <- union(cells, sensitivityInfo(object)[['cellid']])
        }
	}
	drugInfo(object) <- drugInfo(object)[drugs , , drop=drop]
	cellInfo(object) <- cellInfo(object)[cells , , drop=drop]
	object@curation$drug <- object@curation$drug[drugs , , drop=drop]
	object@curation$cell <- object@curation$cell[cells , , drop=drop]
	object@curation$tissue <- object@curation$tissue[cells , , drop=drop]
	if (object@datasetType == 'sensitivity' | object@datasetType == 'both'  & length(exps) == 0) {
	  object@sensitivity$n <- object@sensitivity$n[cells, drugs, drop=drop]
	}
	if (object@datasetType == 'perturbation' | object@datasetType == 'both') {
	  object@perturbation$n <- object@perturbation$n[cells, drugs, , drop=drop]
    }
      return(object)
}


# ==== LongTable Class

## FIXME:: Update this documentation!
#' Subset method for a LongTable object.
#'
#' Allows use of the colData and rowData `data.table` objects to query based on
#'  rowID and colID, which is then used to subset all value data.tables stored
#'  in the dataList slot.
#'
#' This function is endomorphic, it always returns a LongTable object.
#'
#' @param x [`LongTable`] The object to subset.
#' @param i [`character`], [`numeric`], [`logical`] or [`expression`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all row id columns matching the vector.
#'
#'  Numeric or Logical: these select based on the rowKey from the `rowData`
#'      method for the `LongTable`.
#'
#'  Expression: Accepts valid query statements to the `data.table` i parameter,
#'      this can be used to make complex queries using the `data.table` API
#'      for the `rowData` data.table.
#'
#' @param j [`character`], [`numeric`], [`logical`] or [`expression`]
#'  Character: pass in a character vector of drug names, which will subset the
#'      object on all drug id columns matching the vector.
#'
#'  Numeric or Logical: these select based on the rowID from the `rowData`
#'      method for the `LongTable`.
#'
#'  Expression: Accepts valid query statements to the `data.table` i parameter,
#'      this can be used to make complex queries using the `data.table` API
#'      for the `colData` data.table.
#'
#' @param values [`character`, `numeric` or `logical`] Optional list of value
#'      names to subset. Can be used to subset the dataList column further,
#'      returning only the selected items in the new LongTable.
#' @param reindex [`logical`] Should the col/rowKeys be remapped after subsetting.
#'      defaults to FALSE, since reindexing can have significant performance
#'      costs in chained subsetting. We recommend using the `reindex` function
#'      to manually reset the indexes after chained subsets.
#'
#' @return [`LongTable`] A new `LongTable` object subset based on the specified
#'      parameters.
#'
#' @importMethodsFrom BiocGenerics subset
#' @importFrom crayon magenta cyan
#' @import data.table
#' @export
setMethod('subset', signature('LongTable'), function(x, i, j, assays, reindex=FALSE) {

    longTable <- x
    rm(x)

    # local helper functions
    .rowData <- function(...) rowData(..., key=TRUE)
    .colData <- function(...) colData(..., key=TRUE)
    .tryCatchNoWarn <- function(...) suppressWarnings(tryCatch(...))
    .strSplitLength <- function(...) length(strsplit(...))

    # subset rowData
    ## FIXME:: Can I parameterize this into a helper that works for both row
    ## and column data?
    if (!missing(i)) {
        ## TODO:: Clean up this if-else block
        if (.tryCatchNoWarn(is.call(i), error=function(e) FALSE)) {
            # Do nothing
        } else if (.tryCatchNoWarn(is.character(i), error=function(e) FALSE)) {
            ## TODO:: Implement diagnosis for failed regex queries
            idCols <- colnames(.rowIDData(longTable))
            if (max(unlist(lapply(i, .strSplitLength, split=':'))) > length(idCols))
                stop(cyan$bold('Attempting to select more rowID columns than
                    there are in the LongTable.\n\tPlease use query of the form ',
                    paste0(idCols, collapse=':')))
            i <- grepl(.preprocessRegexQuery(i), rownames(longTable), ignore.case=TRUE)
            i <- str2lang(.variableToCodeString(i))
        } else {
            i <- substitute(i)
        }
        rowDataSubset <- .rowData(longTable)[eval(i), ]
    } else {
        rowDataSubset <- .rowData(longTable)
    }

    # subset colData
    if (!missing(j)) {
        ## TODO:: Clean up this if-else block
        if (.tryCatchNoWarn(is.call(j), error=function(e) FALSE, silent=TRUE)) {
            # Do nothing
        } else if (.tryCatchNoWarn(is.character(j), error=function(e) FALSE, silent=TRUE)) {
            ## TODO:: Implement diagnosis for failed regex queries
            idCols <- colnames(.colIDData(longTable))
            if (max(unlist(lapply(j, .strSplitLength, split=':'))) > length(idCols))
                stop(cyan$bold('Attempting to select more ID columns than there
                    are in the LongTable.\n\tPlease use query of the form ',
                    paste0(idCols, collapse=':')))
            j <- grepl(.preprocessRegexQuery(j), colnames(longTable), ignore.case=TRUE)
            j <- str2lang(.variableToCodeString(j))
        } else {
            j <- substitute(j)
        }
        colDataSubset <- .colData(longTable)[eval(j), ]
    } else {
        colDataSubset <- .colData(longTable)
    }

    # Subset assays to only keys in remaining in rowData/colData
    rowKeys <- rowDataSubset$rowKey
    colKeys <- colDataSubset$colKey

    if (missing(assays)) { assays <- assayNames(longTable) }
    keepAssays <- assayNames(longTable) %in% assays

    assayData <- lapply(assays(longTable)[keepAssays],
                     FUN=.filterLongDataTable,
                     indexList=list(rowKeys, colKeys))

    # Subset rowData and colData to only keys contained in remaining assays
    if (!all(assays %in% assayNames(longTable))) {
        ## TODO:: Implement message telling users which rowData and colData
        ## columns are being dropped when selecting a specific assay.
        assayRowIDs <- unique(Reduce(c, lapply(assayData, `$`, name='rowKey')))
        assayColIDs <- unique(Reduce(c, lapply(assayData, `$`, name='colKey')))

        rowDataSubset <- rowDataSubset[rowKey %in% assayRowIDs]
        colDataSubset <- colDataSubset[colKey %in% assayColIDs]
    }

    newLongTable <- LongTable(colData=colDataSubset, colIDs=longTable@.intern$colIDs ,
                     rowData=rowDataSubset, rowIDs=longTable@.intern$rowIDs,
                     assays=assayData, metadata=metadata(longTable))

    return(if (reindex) reindex(newLongTable) else newLongTable)
})

#' Convenience function for converting R code to a call
#'
#' This is used to pass through unevaluated R expressions into subset and
#'   `[`, where they will be evaluated in the correct context.
#'
#' @param ... [`parilist`] Arbitrary R code to subsitute.
#'
#' @return [`call`] An R call object containing the R code from `...`
#'
#' @export
. <- function(...) substitute(...)

# ---- subset LongTable helpers

#' Collapse vector of regex queries with | and replace * with .*
#'
#' @param queryString [`character`] Raw regex queries.
#'
#' @return [`character`] Formatted regex query.
#'
#' @keywords internal
#' @noRd
.preprocessRegexQuery <- function(queryString) {
    # Support vectors of regex queries
    query <- paste0(queryString, collapse='|')
    # Swap all * with .*
    query <- gsub('\\.\\*', '*', query)
    return(gsub('\\*', '.*', query))
}

#'
#'
#'
#'
#' @keywords internal
#' @noRd
.validateRegexQuery <- function(regex, names) {
    ## TODO:: return TRUE if reqex query is valid, otherwise return error message
}

#' Convert an R object in a variable into a string of the code necessary to
#'   create that object
#'
#' @param variable [`Symbol`] A symbol containing an R variable
#'
#' @return [`string`] A string representation of the code necessary to
#'   reconstruct the variable.
#'
#' @keywords internal
#' @noRd
.variableToCodeString <- function(variable) {
    codeString <- paste0(capture.output(dput(variable)), collapse='')
    codeString <- gsub('\"', "'", codeString)
    return(codeString)
}
