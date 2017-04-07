#################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:	
##      - data: gene expression data matrix
##			- drug: single or vector of drug(s) of interest; if a vector of drugs is provided, they will be considered as being the same drug and will be jointly analyszed
##			- drug.id: drug used in each experiment
##			- drug.concentration: drug concentration used in each experiment
##			- type: cell or tissue type for each experiment
##			- xp: type of experiment (perturbation or control)
##      - batch: experiment batches
##      - duration: The duration of the experiment, in a consistent unit
##			- single.type: Should the statitsics be computed for each cell/tissue type separately?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs:
## list of datafraes with the statistics for each gene, for each type
##			- list of data.frame with similar results for each type line separately if any
##
#################################################

rankGeneDrugPerturbation <-
function (data, drug, drug.id, drug.concentration, type, xp, batch, duration, single.type=FALSE, nthread=1, verbose=FALSE) {

	if (nthread != 1) {
		availcore <- parallel::detectCores()
		if (missing(nthread) || nthread < 1 || nthread > availcore) {
			# print(paste("available cores",availcore,"allocated"))
			nthread <- availcore
		}
		else{
			# print(paste("all",nthread,"cores have been allocated"))
		}
	}
	if (any(c(length(drug.id), length(drug.concentration), length(type), length(xp), length(batch), length(duration)) != nrow(data))) {
		stop("length of drug.id, drug.concentration, type, xp, duration and batch should be equal to the number of rows of data!")
	}
	names(drug.id) <- names(drug.concentration) <- names(type) <- names(batch) <- names(duration) <- rownames(data)
	if (!all(complete.cases(type, xp, batch, duration))) {
		stop("type, batch, duration and xp should not contain missing values!")
	}
## is the drug in the dataset?
	drugix <- drug.id %in% drug
		  
	if (sum(drugix) == 0) {
		warning(sprintf("Drug(s) %s not in the dataset", paste(drug, collapse=", ")))
		return(list("all.type"=NULL, "single.type"=NULL))
	}
## select xps with controls or with the drug(s) of interest
	iix <- xp=="control" | drugix
	data <- data[iix, ,drop=FALSE]
	drug.id <- drug.id[iix]
	drug.concentration <- drug.concentration[iix]
	type <- type[iix]
	xp <- xp[iix]
	batch <- batch[iix]
	duration <- duration[iix]
	
	res.type <- NULL
	
## build input matrix
	inpumat <- NULL
## for each batch/vehicle of perturbations+controls (test within each batch/vehicle to avoid batch effect)
	ubatch <- sort(unique(batch[!is.na(xp) & xp == "perturbation"]))
	names(ubatch) <- paste("batch", ubatch, sep="")
	
	for (bb in 1:length(ubatch)) {
## identify the perturbations and corresponding control experiments
		xpix <- rownames(data)[complete.cases(batch, xp) & batch == ubatch[bb] & xp == "perturbation"]
		ctrlix <- rownames(data)[complete.cases(batch, xp) & batch == ubatch[bb] & xp == "control"]
		
		if (all(!is.na(c(xpix, ctrlix))) && length(xpix) > 0 && length(ctrlix) > 0) {
			if (!all(is.element(ctrlix, rownames(data)))) {
				stop("data for some control experiments are missing!")
			}
			if (verbose) {
				cat(sprintf("type %s: batch %i/%i -> %i vs %i\n", utype[bb], bb, length(ubatch), length(xpix), length(ctrlix)))
			}
## transformation of drug concentrations values
			conc <- drug.concentration * 10^6
			inpumat <- rbind(inpumat, data.frame("treated"=c(rep(1, length(xpix)), rep(0, length(ctrlix))), "type"=c(type[xpix], type[ctrlix]), "batch"=paste("batch", c(batch[xpix], batch[ctrlix]), sep=""), "concentration"=c(conc[xpix], conc[ctrlix]), "duration"= c(duration[xpix], duration[ctrlix])))
		}
	}

	inpumat[ , "type"] <- factor(inpumat[ , "type"], ordered=FALSE)
	inpumat[ , "batch"] <- factor(inpumat[ , "batch"], ordered=FALSE)

	if (nrow(inpumat) < 3 || length(sort(unique(inpumat[ , "concentration"]))) < 2 || length(unique(inpumat[ , "duration"])) < 2) {
## not enough experiments in drug list
		warning(sprintf("Not enough data for drug(s) %s", paste(drug, collapse=", ")))
		return(list("all.type"=NULL, "single.type"=NULL))
	}
	
	res <- NULL
	utype <- sort(unique(as.character(inpumat[ , "type"])))
	ltype <- list("all"=utype)
	if(single.type) {
		ltype <- c(ltype, as.list(utype))
		names(ltype)[-1] <- utype
	}
	for(ll in 1:length(ltype)) {
## select the type of cell line/tissue of interest
		inpumat2 <- inpumat[!is.na(inpumat[ , "type"]) & is.element(inpumat[ , "type"], ltype[[ll]]), , drop=FALSE]
		inpumat2 <- inpumat2[complete.cases(inpumat2), , drop=FALSE]
		if (nrow(inpumat2) < 3 || length(sort(unique(inpumat2[ , "concentration"]))) < 2) {
## not enough experiments in data
			nc <- c("estimate", "se", "n", "tstat", "fstat", "pvalue")
			rest <- matrix(NA, nrow=nrow(data), ncol=length(nc), dimnames=list(rownames(data), nc))
		} else {
## test perturbation vs control
			if(nthread > 1) {
## parallel threads
				splitix <- parallel::splitIndices(nx=ncol(data), ncl=nthread)
				splitix <- splitix[sapply(splitix, length) > 0]
				mcres <- parallel::mclapply(splitix, function(x, data, inpumat) {
											res <- t(apply(data[rownames(inpumat), x, drop=FALSE], 2, geneDrugPerturbation, concentration=inpumat[ , "concentration"], type=inpumat[ , "type"], batch=inpumat[ , "batch"], duration=inpumat[,"duration"]))
											return(res)
											}, data=data, inpumat=inpumat2)
				rest <- do.call(rbind, mcres)
			} else {
				rest <- t(apply(data[rownames(inpumat2), , drop=FALSE], 2, geneDrugPerturbation, concentration=inpumat2[ , "concentration"], type=inpumat2[ , "type"], batch=inpumat2[ , "batch"], duration=inpumat2[,"duration"]))
			}
		}
		rest <- cbind(rest, "fdr"=p.adjust(rest[ , "pvalue"], method="fdr"))
		res <- c(res, list(rest))
	}
	names(res) <- names(ltype)
	return(res)
}

## End
