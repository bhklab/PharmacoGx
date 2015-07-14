########################
## Petr Smirnov
## All rights Reserved
## May 17, 2015
## Function to intersect objects of the PharmacoSet class type
########################


#' Intersects objects of the PharmacoSet class, subsetting them to the common genes, drugs and/or cell lines as selected by the user.
#' 
#' Given a list of PharmacoSets, the function will find the common genes, drugs,
#' and/or cell lines, and return PharmacoSets that contain data only pertaining 
#' to the common genes, drugs, and/or cell lines. The mapping between dataset 
#' drug, cell and gene names is done either using annotations found in the 
#' PharmacoSet object's internal curation slot, or according to matching tables 
#' provided by the user
#' 
#' @examples 
#' data(CGPsmall)
#' data(CCLEsmall)
#' common <- intersectPSet(list('CGP'=CGPsmall, 'CCLE'=CCLEsmall),
#'                         intersectOn = c("drugs", "cell.lines"))
#' common$CGP
#' common$CCLE
#' 
#' 
#' @param pSets [list] a list of PharmacoSet objects, of which the function
#'   should find the intersection
#' @param intersectOn [character] which identifiers to intersect on, genes,
#'   drugs or cell lines
#' @param drugMatch [matrix] a matrix which specifies the mapping between each 
#'   dataset in the list. Should have the same number of columns as objects in 
#'   the pSets list. Each column should have the Drug Idetifiers used in that
#'   pSet. Each row will then map the different identifiers across datasets.
#'   The name of each row should be set to the universal identifier you wish
#'   to use across all the datasets. If no drugMatch table is provided, the
#'   function will attempt to use the curation slot of each pSet to create 
#'   a mapping between the datasets. 
#' @param cellMatch [matrix] a matrix which specifies the mapping between each 
#'   dataset in the list. Should have the same number of columns as objects in 
#'   the pSets list. Each column should have the Cell Idetifiers used in that
#'   pSet. Each row will then map the different identifiers across datasets.
#'   The name of each row should be set to the universal identifier you wish
#'   to use across all the datasets. If no cellMatch table is provided, the
#'   function will attempt to use the curation slot of each pSet to create 
#'   a mapping between the datasets. 
#' @param geneMatch [matrix] a matrix which specifies the mapping between each 
#'   dataset in the list. Should have the same number of columns as objects in 
#'   the pSets list. Each column should have the Gene/Probe Idetifiers used in that
#'   pSet. Each row will then map the different identifiers across datasets.
#'   The name of each row should be set to the universal identifier you wish
#'   to use across all the datasets. If no geneMatch table is provided, the
#'   function will assume that the gene/probe names already match between datasets.
#' @param strictIntersect [boolean] Should the intersection keep only the drugs 
#'   and cell lines that have been tested on together?
#' @param verbose [boolean] Should the function announce its key steps?
#' @return [list] a list of pSets, contatining only the intersection
#' @export
#' 
intersectPSet <- function (
    pSets=list(),
    intersectOn=c("drugs","genes","cell.lines", "concentrations"), 
    drugMatch=NULL, 
    cellMatch=NULL,
    geneMatch=NULL,
    strictIntersect=FALSE,
    verbose=TRUE) {
  
  ## TODO: Fix the strict intersection!!!!!!
	if (length(pSets) < 1){ stop("No PharmacoSets passed to function")}
	if (length(pSets) == 1){
		pSet <- pSets[[1]]
		# cellx <- match(cellMatch[,1], rownames(cellInfo(pSet)))
		# drugx <- match(drugMatch[,1], rownames(drugInfo(pSet)))
		# genex <- match(geneMatch[,1], rownames(geneInfo(pSet)))
		if (all(c("drugs", "cell.lines", "genes") %in% intersectOn)){
			### Note, this is NOT redundant, as the current implementation guarantees that the drugs and cells passed to subset will be in the drug 
			### and cell tables even if they do not exist anymore in dataset. This is so we can update the names correctly in the data sets.
			pSet <- subsetTo(pSet, drugs=drugMatch[,1], cells=cellMatch[,1], genes=geneMatch[,1])
			pSet <- updateDrugId(pSet, rownames(drugMatch))
			pSet <- updateCellId(pSet, rownames(cellMatch))
		} else {
			if ((all(c("drugs", "cell.lines")%in% intersectOn ))){
			### Note, this is NOT redundant, as the current implementation guarantees that the drugs and cells passed to subset will be in the drug 
			### and cell tables even if they do not exist anymore in dataset. This is so we can update the names correctly in the data sets.
				pSet <- subsetTo(pSet, drugs=drugMatch[,1], cells=cellMatch[,1])
				pSet <- updateDrugId(pSet, rownames(drugMatch))
				pSet <- updateCellId(pSet, rownames(cellMatch))
			} else {
				if ("drugs" %in% intersectOn) {
					pSet <- subsetTo(pSet, drugs=drugMatch[,1])
					pSet <- updateDrugId(pSet, rownames(drugMatch))
				} 
				if ("cell.lines"%in% intersectOn){
					pSet <- subsetTo(pSet, cells=cellMatch[,1])
					pSet <- updateCellId(pSet, rownames(cellMatch))
				} 
				if ("genes"%in% intersectOn){
					pSet <- subsetTo(pSet, genes=geneMatch[,1])
					rownames(geneInfo(pSet)) <- rownames(geneMatch)
				}
			}
		}
		return(pSet)
	}
	if (length(pSets) > 1){

    temp.names <- names(pSets)
		if ("drugs" %in% intersectOn){

			if(is.null(drugMatch)){
				# common.drugs <- intersectList(lapply(pSets, drugNames))
				# if(length(common.drugs)==0){stop("No common drugs in all datasets")}
				# drugMatch <- matrix(rep(common.drugs, times=length(pSets)),nrow=length(common.drugs), ncol=length(pSets))
        		# rownames(drugMatch) <- common.drugs
				common.drugs <- intersectList(lapply(pSets, function(x){return(x@curation$drug[,"unique.drugid"])}))
				drugMatch <- data.frame(lapply(pSets, function(x, common.drugs){return(rownames(x@curation$drug)[match(common.drugs,x@curation$drug[,"unique.drugid"])])}, common.drugs=common.drugs))
				drugMatch <- as.matrix(drugMatch)
				rownames(drugMatch) <- common.drugs
			}
		} 
		if ("cell.lines" %in% intersectOn){

			if(is.null(cellMatch)){
				# common.cells <- intersectList(lapply(pSets, cellNames))
				# if(length(common.cells)==0){stop("No common cell lines in all datasets")}
				# cellMatch <- matrix(rep(common.cells, times=length(pSets)),nrow=length(common.cells), ncol=length(pSets))
        		# rownames(cellMatch) <- common.cells
				common.cells <- intersectList(lapply(pSets, function(x){return(x@curation$cell[,"unique.cellid"])}))
				cellMatch <- data.frame(lapply(pSets, function(x, common.cells){return(rownames(x@curation$cell)[match(common.cells,x@curation$cell[,"unique.cellid"])])}, common.cells=common.cells))
				cellMatch <- as.matrix(cellMatch)
				rownames(cellMatch) <- common.cells
			}
		} 

		if ("genes" %in% intersectOn){

			if(is.null(geneMatch)){
				common.genes <- intersectList(lapply(pSets, function(x,...){return(rownames(geneInfo(x)))}))
				if(length(common.genes)==0){stop("No common genes in all datasets")}
				geneMatch <- matrix(rep(common.genes, times=length(pSets)),nrow=length(common.genes), ncol=length(pSets))
        rownames(geneMatch) <- common.genes

			}
		}
		
		pSets <- lapply(1:length(pSets), function(x, pSets, intersectOn, drugMatch, cellMatch, geneMatch, verbose){
    
      return(intersectPSet(pSets[x], intersectOn, drugMatch[,x,drop=FALSE], cellMatch[,x,drop=FALSE], geneMatch[,x,drop=FALSE], verbose=verbose))
    }, pSets=pSets, intersectOn=intersectOn, drugMatch=drugMatch, cellMatch=cellMatch, geneMatch=geneMatch, verbose=verbose)
	
	# if ("concentrations"%in% intersectOn){
#         pSets <- .intesectConcentrations(pSets,cap=100)
#     }
   
      names(pSets) <-  temp.names     
      if ("concentrations" %in% intersectOn){
        pSets <- .calculateSensitivities(pSets, cellMatch, drugMatch, cap=100)
      }

  return(pSets)
	}
}


# .intesectConcentrations <- function(pSets=list(), cap=100){
	
#     # pSet <- .calculateSensitivities(pSet, dose.range=dose.range, cap=cap)
	
	
# 	cellnall <- intersect(lapply(pSets, cellNames))
# 	drugnall <- intersect(lapply(pSets, drugNames))
	
# 	exp_row <- max(lapply(pSets, function(x){
# 		return(nrow(x@sensitivity$info))
# 	}))
	
# 	expmatch <- matrix(NA, nrow=exp_row, ncol=length(pSets))
# 	row = 1
# 	for (drug in drugnall){		
# 		for (cell in cellnall){
		
# 			exps <- lapply(pSets, function(x){
				
# 				sI <- sensitivityInfo(x)
				
# 				return(rownames(sI[which(sI[,"drugid"]==drug & sI[,"cellid"]==cell),]))
				
# 			})
# 			if (max(length(exps))==1){
# 				exp_match[i,] <- unlist(exps)
# 			}
# 			### TODO:: deal with duplicates
# 		}
			
# 		}
		
	
	
	
# 	### TODO:: use experiment matching columns to find common concentration ranges and match to that. How do we ensure that for replicates we dont end up
# 	### reducing concetration range if the concentration ranges of the replicates are unequal??? How do we deal with that???
	
	
# 	min_conc <- pmin(lapply(pSets, function(x){
# 		return(sensitivityInfo(x)[,"min.Dose.nM"])
# 	}))
# 	max_conc <- pmax(lapply(pSets, function(x){
# 		return(sensitivityInfo(x)[,"max.Dose.nM"])
# 	}))
	

	
# }
