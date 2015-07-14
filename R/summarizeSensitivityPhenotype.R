########################
## Petr Smirnov
## All rights Reserved
## June 12, 2015
########################


#' Takes the sensitivity data from a PharmacoSet, and summarises them into a
#' drug vs cell line table
#' 
#' This function creates a table with cell lines as rows and drugs as columns,
#' summarising the drug senstitivity data of a PharmacoSet into drug-cell line
#' pairs
#' 
#' @examples 
#' data(CGPsmall)
#' CGPauc <- summarizeSensitivityPhenotype(CGPsmall, sensitivity.measure='auc_published')
#'
#' @param pSet [PharmacoSet] The PharmacoSet from which to extract the data
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? The current
#'   choices are 'ic50_published', 'auc_published', 'ic50_recomputed',
#'   'auc_recomputed'.
#' @param summaryStat [character] which summaryStat method to use if there are repeated
#'   pairs of cells and drugs?
#' @return [matrix] A matrix with cell lines going down the rows, drugs across
#'   the columns, with the selected sensitivity statistic for each pair.
#' @export


summarizeSensitivityPhenotype <- function(pSet, sensitivity.measure=c("ic50_published", "auc_published", "ic50_recompted", "auc_recomputed"), summaryStat=c("median","mean", "first", "last")){
	
	summaryStat <- match.arg(summaryStat)
  sensitivity.measure <- match.arg(sensitivity.measure)
  
	cellnall <- cellNames(pSet)
	drugnall <- drugNames(pSet)
	
	dd <- matrix(NA, nrow=length(cellnall), ncol=length(drugnall), dimnames=list(cellnall, drugnall))
	
	sensInfo <- sensitivityInfo(pSet) 
	sensData <- sensitivityData(pSet)
	
	for (i in 1:length(drugnall)){
		drugX <- sensInfo[,"drugid"]==drugnall[i]
		
		for (j in 1:length(cellnall)){
			rowX <- drugX & (sensInfo[,"cellid"]==cellnall[j])
			if (sum(rowX)==1){
				dd[j,i] <- sensData[rowX, sensitivity.measure]
			}
			if (sum(rowX)>1){
				switch(summaryStat, "mean"={
					
					dd[j,i] <- mean(sensData[rowX, sensitivity.measure])
					
				},
				"median"={
					
					dd[j,i] <- median(sensData[rowX, sensitivity.measure])
					
				}, 
				"first"={
					
					dd[j,i] <- sensData[rowX[1], sensitivity.measure]
					
				},
				"last" = {
					
					dd[j,i] <- sensData[rowX[length(rowX)], sensitivity.measure]
					
				})
			}
			
		}
		
	}
	
  return(dd)
	
}



# setGeneric("updateTables", function(pSet, summaryStat, sensitivity.measures, verbose) standardGeneric("updateTables"))
# setMethod(updateTables, "PharmacoSet", function(pSet, summaryStat=vector(), sensitivity.measures=vector(), verbose=FALSE){
# 	# summaryStat=c("median","mean", "first", "last")
# 	# summaryStat=match.arg(summaryStat)
# 	if (length(sensitivity.measures)==0){
# 		sensitivity.measures=sensitivitysensitivity.measures(pSet)
# 	}
# 	if (length(summaryStat)==0){
# 		summaryStat = rep("median", times=length(sensitivity.measures))
# 	} else if (length(summaryStat)!=length(sensitivity.measures)){
# 		stop("The length of the summaryStat argument and the sensitivity.measures argument do not match")
# 	}
#
#     for (i in length(sensitivity.measures)){
#     	  if(verbose){
#   		  print("Now generating table for: ")
#   		  print(sensitivity.measures[[i]])
#     	  }
#   	  pSet@tables[,,sensitivity.measures[[i]]] <- summarizeSensitivityPhenotype(pSet, sensitivity.measure=sensitivity.measures[[i]], summaryStat=summaryStat[[i]])
#
#     }
#
# 	pSet@table.summaryStat[[i]] <- summaryStat[[i]]
#
# 	return(pSet)
#
# })
#
# setGeneric("sensitivityTable", function(pSet, sensitivity.measure) standardGeneric("sensitivityTable"))
# setMethod(sensitivityTable, "PharmacoSet", function(pSet, sensitivity.measure="auc_published"){
#
# 	return(matrix(pSet@tables[,,sensitivity.measure, drop=FALSE], nrow=dim(pSet@tables)[[1]], ncol=dim(pSet@tables)[[2]], dimnames=dimnames(pSet@tables)[1:2]))
#
# })