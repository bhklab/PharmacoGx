#' @include PharmacoSet-class.R
NULL

mergePSets <- function(mDataPSet, sensDataPSet, commonCellsOnly=FALSE, ...){

	if(commonCellsOnly){
		common <- intersectPSet(list(mDataPSet, sensDataPSet), intersectOn=c("celllines"))
		mDataPSet <- common[[1]]
		sensDataPSet <- common[[2]]
	}

  ### union of cell lines
	ucell <- union(sampleNames(sensDataPSet), sampleNames(mDataPSet))

  ### molecular profiles
	mergePSet <- sensDataPSet
	molecularProfilesSlot(mergePSet) <- molecularProfilesSlot(mDataPSet)

  ### number of sensitivity experiments
  #acell <- setdiff(rownames(sampleInfo(mDataPSet)), rownames(sampleInfo(sensDataPSet)))

	### cell line annotations
	new.cell.info <- c(union(colnames(sampleInfo(sensDataPSet)), colnames(sampleInfo(mDataPSet))), "dataset")
	cell.info.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.info), dimnames=list(ucell, new.cell.info)), check.names=FALSE)
	cell.info.df[rownames(sampleInfo(mDataPSet)), colnames(sampleInfo(mDataPSet))] <- sampleInfo(mDataPSet)
	cell.info.df[rownames(sampleInfo(sensDataPSet)), colnames(sampleInfo(sensDataPSet))] <- sampleInfo(sensDataPSet)
	cell.info.df[setdiff(rownames(sampleInfo(sensDataPSet)), rownames(sampleInfo(mDataPSet))), "dataset"] <- name(sensDataPSet)
	cell.info.df[setdiff(rownames(sampleInfo(mDataPSet)), rownames(sampleInfo(sensDataPSet))), "dataset"] <- name(mDataPSet)
	cell.info.df[intersect(rownames(sampleInfo(mDataPSet)), rownames(sampleInfo(sensDataPSet))), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	sampleInfo(mergePSet) <- cell.info.df

	### curation of cell line names
	new.cell.curation <- c(union(colnames(curation(sensDataPSet)$sample), colnames(curation(mDataPSet)$sample)), "dataset")
	cell.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.cell.curation), dimnames=list(ucell, new.cell.curation)), check.names=FALSE)
	cell.curation.df[rownames(curation(mDataPSet)$sample), colnames(curation(mDataPSet)$sample)] <- curation(mDataPSet)$sample
	cell.curation.df[rownames(curation(sensDataPSet)$sample), colnames(curation(sensDataPSet)$sample)] <- curation(sensDataPSet)$sample
	cell.curation.df[setdiff(rownames(curation(sensDataPSet)$sample), rownames(curation(mDataPSet)$sample)), "dataset"] <- name(sensDataPSet)
	cell.curation.df[setdiff(rownames(curation(mDataPSet)$sample), rownames(curation(sensDataPSet)$sample)), "dataset"] <- name(mDataPSet)
	cell.curation.df[intersect(rownames(curation(mDataPSet)$sample), rownames(curation(sensDataPSet)$sample)), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	curation(mergePSet)$sample <- cell.curation.df

	### curation of tissue names
	new.tissue.curation <- c(union(colnames(curation(sensDataPSet)$tissue), colnames(curation(mDataPSet)$tissue)), "dataset")
	tissue.curation.df <- data.frame(matrix(NA, nrow=length(ucell), ncol=length(new.tissue.curation), dimnames=list(ucell, new.tissue.curation)), check.names=FALSE)
	tissue.curation.df[rownames(curation(mDataPSet)$tissue), colnames(curation(mDataPSet)$tissue)] <- curation(mDataPSet)$tissue
	tissue.curation.df[rownames(curation(sensDataPSet)$tissue), colnames(curation(sensDataPSet)$tissue)] <- curation(sensDataPSet)$tissue
	tissue.curation.df[setdiff(rownames(curation(sensDataPSet)$tissue), rownames(curation(mDataPSet)$tissue)), "dataset"] <- name(sensDataPSet)
	tissue.curation.df[setdiff(rownames(curation(mDataPSet)$tissue), rownames(curation(sensDataPSet)$tissue)), "dataset"] <- name(mDataPSet)
	tissue.curation.df[intersect(rownames(curation(mDataPSet)$tissue), rownames(curation(sensDataPSet)$tissue)), "dataset"] <- paste0(name(sensDataPSet), "///", name(mDataPSet))
	curation(mergePSet)$tissue <- tissue.curation.df

	sensNumber(mergePSet) <- PharmacoGx:::.summarizeSensitivityNumbers(mergePSet)

	annotation(mergePSet)$name <- paste(name(mDataPSet), name(sensDataPSet), sep=".")

	annotation(mergePSet)$dateCreated <- date()

	checkPSetStructure(mergePSet)

    return(mergePSet)
}
