

mergePSets <- function(mDataPSet, sensDataPSet, commonCellsOnly=FALSE, ...){

	if(commonCellsOnly){
		common <- intersectPSet(list(mDataPSet, sensDataPSet), intersectOn=c("celllines"))
		mDataPSet <- common[[1]]
		sensDataPSet <- common[[2]]
	}

  ### union of cell lines
	ucell <- union(cellNames(sensDataPSet), cellNames(mDataPSet))

  ### molecular profiles
	mergePSet <- sensDataPSet
	mergePSet@molecularProfiles <- mDataPSet@molecularProfiles

  ### number of sensitivity experiments
  #acell <- setdiff(rownames(cellInfo(mDataPSet)), rownames(cellInfo(sensDataPSet)))

	### cell line annotations
	new.cell.info <- c(union(colnames(cellInfo(sensDataPSet)), colnames(cellInfo(mDataPSet))), "dataset")
	cell.info.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.info), dimnames=list(ucell, new.cell.info)), check.names=FALSE)
	cell.info.df[rownames(cellInfo(mDataPSet)), colnames(cellInfo(mDataPSet))] <- cellInfo(mDataPSet)
	cell.info.df[rownames(cellInfo(sensDataPSet)), colnames(cellInfo(sensDataPSet))] <- cellInfo(sensDataPSet)
	cell.info.df[setdiff(rownames(cellInfo(sensDataPSet)), rownames(cellInfo(mDataPSet))), "dataset"] <- name(sensDataPSet)
	cell.info.df[setdiff(rownames(cellInfo(mDataPSet)), rownames(cellInfo(sensDataPSet))), "dataset"] <- name(mDataPSet)
	cell.info.df[intersect(rownames(cellInfo(mDataPSet)), rownames(cellInfo(sensDataPSet))), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	cellInfo(mergePSet) <- cell.info.df

	### curation of cell line names
	new.cell.curation <- c(union(colnames(sensDataPSet@curation$cell), colnames(mDataPSet@curation$cell)), "dataset")
	cell.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.curation), dimnames=list(ucell, new.cell.curation)), check.names=FALSE)
	cell.curation.df[rownames(mDataPSet@curation$cell), colnames(mDataPSet@curation$cell)] <- mDataPSet@curation$cell
	cell.curation.df[rownames(sensDataPSet@curation$cell), colnames(sensDataPSet@curation$cell)] <- sensDataPSet@curation$cell
	cell.curation.df[setdiff(rownames(sensDataPSet@curation$cell), rownames(mDataPSet@curation$cell)), "dataset"] <- name(sensDataPSet)
	cell.curation.df[setdiff(rownames(mDataPSet@curation$cell), rownames(sensDataPSet@curation$cell)), "dataset"] <- name(mDataPSet)
	cell.curation.df[intersect(rownames(mDataPSet@curation$cell), rownames(sensDataPSet@curation$cell)), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	mergePSet@curation$cell <- cell.curation.df

	### curation of tissue names
	new.tissue.curation <- c(union(colnames(sensDataPSet@curation$tissue), colnames(mDataPSet@curation$tissue)), "dataset")
	tissue.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.tissue.curation), dimnames=list(ucell, new.tissue.curation)), check.names=FALSE)
	tissue.curation.df[rownames(mDataPSet@curation$tissue), colnames(mDataPSet@curation$tissue)] <- mDataPSet@curation$tissue
	tissue.curation.df[rownames(sensDataPSet@curation$tissue), colnames(sensDataPSet@curation$tissue)] <- sensDataPSet@curation$tissue
	tissue.curation.df[setdiff(rownames(sensDataPSet@curation$tissue), rownames(mDataPSet@curation$tissue)), "dataset"] <- name(sensDataPSet)
	tissue.curation.df[setdiff(rownames(mDataPSet@curation$tissue), rownames(sensDataPSet@curation$tissue)), "dataset"] <- name(mDataPSet)
	tissue.curation.df[intersect(rownames(mDataPSet@curation$tissue), rownames(sensDataPSet@curation$tissue)), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	mergePSet@curation$tissue <- tissue.curation.df

	mergePSet@treatmentResponse$n <- PharmacoGx:::.summarizeSensitivityNumbers(mergePSet)

	mergePSet@annotation$name <- paste(name(mDataPSet), name(sensDataPSet), sep=".")

	mergePSet@annotation$dateCreated <- date()

	checkPSetStructure(mergePSet)

    return(mergePSet)
}
