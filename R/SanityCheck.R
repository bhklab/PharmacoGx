sanitizeInput <- function(conc,
	viability,
	Hill_fit,
	conc_as_log = FALSE,
	viability_as_pct = TRUE,
	trunc = TRUE,
	verbose = TRUE # Set to 2 to see debug printouts
	){


	if (is.logical(conc_as_log) == FALSE) {
		print(conc_as_log)
		stop("'conc_as_log' is not a logical.")
	}

	if (is.logical(viability_as_pct) == FALSE) {
		print(viability_as_pct)
		stop("'viability_as_pct' is not a logical.")
	}

	if (is.logical(trunc) == FALSE) {
		print(trunc)
		stop("'trunc' is not a logical.")
	}
	if(!is.finite(verbose)){
		stop("'verbose' should be a logical (or numerical) argument.")
	}
	if(!missing(viability)&&!missing(conc)&&missing(Hill_fit))
	{ 
	  if (length(conc) != length(viability)) {
	    if(verbose==2){
	      print(conc)
	      print(viability) 
	    }
	    stop("Log concentration vector is not of same length as viability vector.")
	  }
		if( any(is.na(conc)&(!is.na(viability)))){
			warning("Missing concentrations with non-missing viability values encountered. Removing viability values correspoding to those concentrations")

			myx <- !is.na(conc)
			conc <- as.numeric(conc[myx])
			viability <- as.numeric(viability[myx])

		} 
		if(any((!is.na(conc))&is.na(viability))){

			warning("Missing viability with non-missing concentrations values encountered. Removing concentrations values correspoding to those viabilities")

			myx <- !is.na(viability)
			conc <- as.numeric(conc[myx])
			viability <- as.numeric(viability[myx])

		}
		conc <- as.numeric(conc[!is.na(conc)])
		viability <- as.numeric(viability[!is.na(viability)])
		

  #CHECK THAT FUNCTION INPUTS ARE APPROPRIATE
		if (prod(is.finite(conc)) != 1) {
			print(conc)
			stop("Concentration vector contains elements which are not real numbers.")
		}

		if (prod(is.finite(viability)) != 1) {
			print(viability)
			stop("Viability vector contains elements which are not real numbers.")
		}
		

		if (min(viability) < 0) {
			if (verbose) {
				warning("Warning: Negative viability data.")
			}
		}

		if (max(viability) > (1 + 99 * viability_as_pct)) {
			if (verbose) {
				warning("Warning: Viability data exceeds negative control.")
			}
		}


		if (conc_as_log == FALSE && min(conc) < 0) {
			if (verbose == 2) {
				print(conc)
				print(conc_as_log)
			}
			stop("Negative concentrations encountered. Concentration data may be inappropriate, or 'conc_as_log' flag may be set incorrectly.")
		}

		if (viability_as_pct == TRUE && max(viability) < 5) {
			warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
			if (verbose == 2) {
				print(viability)
				print(viability_as_pct)
			}
		}

		if (viability_as_pct == FALSE && max(viability) > 5) {
			warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
			if (verbose == 2) {
				print(viability)
				print(viability_as_pct)	
			}
		}

		if(is.unsorted(conc)){
			warning("Concentration Values were unsorted. Sorting concentration and ordering viability in same order")
			myx <- order(conc)
			conc <- conc[myx]
			viability <- viability[myx]
		}

  #CONVERT DOSE-RESPONSE DATA TO APPROPRIATE INTERNAL REPRESENTATION
		if (conc_as_log == FALSE ) {
		  ii <- which(conc == 0)
		  if(length(ii) > 0) {
		    conc <- conc[-ii]
		    viability <- viability[-ii]
		  }
			log_conc <- log10(conc)
		} else {
			log_conc <- conc
		}

		if (viability_as_pct == TRUE) {
			viability <- viability / 100
		}
		if (trunc) {
			viability = pmin(as.numeric(viability), 1)
			viability = pmax(as.numeric(viability), 0)
		}

		return(list("log_conc"=log_conc, "viability"=viability))
	} 
	if(!missing(Hill_fit) && missing(viability)){
		if(is.list(Hill_fit)){

			Hill_fit <- unlist(Hill_fit)
		}
		if (conc_as_log == FALSE && Hill_fit[[3]] < 0) {
			print("EC50 passed in as:")
			print(Hill_fit[[3]])
			stop("'conc_as_log' flag may be set incorrectly, as the EC50 is negative when positive value is expected.")
		}

		
		if (viability_as_pct == FALSE && Hill_fit[[2]] > 1) {
			print("Einf passed in as:")
			print(Hill_fit[[2]])
			
			warning("Warning: 'viability_as_pct' flag may be set incorrectly.")
			
		}
		if (conc_as_log == FALSE){
			Hill_fit[[3]] <- log10(Hill_fit[[3]])
		}
		if (viability_as_pct == TRUE){
			Hill_fit[[2]] <- Hill_fit[[2]]/100
		}
		if(missing(conc)){
			return(list("Hill_fit"=Hill_fit))
		} else {
			conc <- as.numeric(conc[!is.na(conc)])
			
			if (prod(is.finite(conc)) != 1) {
				print(conc)
				stop("Concentration vector contains elements which are not real numbers.")
			}
			if (conc_as_log == FALSE && min(conc) < 0) {
				print(conc)
				print(conc_as_log)
				stop("Negative concentrations encountered. Concentration data may be inappropriate, or 'conc_as_log' flag may be set incorrectly.")
			}

			if (conc_as_log == FALSE ) {
				ii <- which(conc == 0)
				if(length(ii) > 0) {
					conc <- conc[-ii]
				}
				log_conc <- log10(conc)
			} else {
				log_conc <- conc
			}
			if(is.unsorted(conc)){
				myx <- order(conc)
				conc <- conc[myx]
			}
			return(list("Hill_fit"=Hill_fit, "log_conc" = log_conc))
		}
		

	}
	if(!missing(Hill_fit)&&!missing(viability)){

		stop("Please pass in only one of 'Hill_fit' and 'viability', it is unclear which to use in the computation.")
	}
	if(missing(Hill_fit)&&missing(viability)){

		stop("Both 'Hill_fit' and 'viability' missing, please pass in some data!")
	}
}
