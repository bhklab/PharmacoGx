
setGeneric("plot")

#' Plots a PharmacoSig object into a Volcano Plot
#' 
#' Given a PharmacoSig, this will plot a volcano plot, with parameters to set cutoffs 
#' for a significant effect size, p value, to pick multiple testing correction strategy, 
#' and to change point colors. Built on top of ggplot, it will return the plot object which
#' can be easily customized as any other ggplot. 
#' 
#' @examples
#' data(GDSCsmall)
#' drug.sensitivity <- drugSensitivitySig(GDSCsmall, mDataType="rna", 
#'              nthread=1, features = fNames(GDSCsmall, "rna")[1])
#' plot(drug.sensitivity)
#' 
#' @param x [PharmacoSig] a PharmacoSig object, result of drugSensitivitySig
#'  or drugPerturbationSig
#' @param adjust.method [character] or [boolean] either FALSE for no adjustment,
#' or one of the methods implemented by p.adjust. Defaults to FALSE for no 
#' correction
#' @param drugs [character] a vector of drug names for which to plot the estimated
#' associations with gene expression 
#' @param features [character] a vector of features for which to plot the estimated
#' associations with drug treatment 
#' @param effect_cutoff the cutoff to use for coloring significant effect sizes. 
#' @param signif_cutoff the cutoff to use for coloring significance by p value or
#' adjusted p values. Not on log scale.
#' @param color one color if no cutoffs set for plotting. A vector of colors otherwise
#' used to color points the in three categories above. 
#' @param ... additional arguments, not currently used, but left here for consistency with plot
#' @return returns a ggplot object, which by default will be evaluated and the plot displayed, or
#' can be saved to a variable for further customization by adding ggplot elements to the returned
#' graph
#' @export
#' @import ggplot2
#' @include signatureClass.R
#' @method plot PharmacoSig
plot.PharmacoSig <- function(x, adjust.method, drugs, features, effect_cutoff, signif_cutoff, color, ...){
	dots <- list(...)
	ndots <- length(dots)
	
	# if(length(dim(x))==2){
	# 	dim(x) <- c(1, dim(x))
	# } else if(length(dim(x)) == 1) {
	# 	dim(x) <- c(1, 1, dim(x))
	# }

	if(missing(adjust.method)){
		adjust.method <- FALSE
	}

	if(missing(drugs)){
		drugs <- colnames(x)
	}

	if(missing(features)){
		features <- rownames(x)
	}

	if(!missing(color)){
		if(!is.null(dots[["colour"]])){
			warning("Both color and colour parameters provided. Will take union of both. This is probably a mistake.")
			color <- union(color, dots[["colour"]])
		}
	} else if (!is.null(dots[["colour"]])){
		color <- dots[["colour"]]
	} ## Case if both missing handled in logic below


	if(isFALSE(adjust.method)){
		p.adjust.f <- function(x) return(x)
	} else {
		p.adjust.f <- function(x) return(p.adjust(x, method=adjust.method))
	}

	x.m <- data.frame(X = as.vector(x[features,drugs,c("estimate")]), 
					  Y = -log10(p.adjust.f(as.vector(x[features,drugs,c("pvalue")]))))
	axis.labs <- c("Estimate", ifelse(isFALSE(adjust.method) || adjust.method == "none", "-Log10 P Value", "-Log10 Corrected P Value"))


	plot.elements <- ggplot() + xlab(axis.labs[1]) + ylab(axis.labs[2])

	
	if(!missing(effect_cutoff) | !missing(signif_cutoff)) {
		x.m$Cutoff <- "Not Significant"

		if(!missing(signif_cutoff)){
			x.m$Cutoff[x.m$Y >= -log10(signif_cutoff)] <- "Significant P Value"

			if(!missing(effect_cutoff)){
				x.m$Cutoff[(x.m$Y >= -log10(signif_cutoff)) & (abs(x.m$X) >= effect_cutoff)] <- "Significant P Value and Effect"
			}

		} else {
			x.m$Cutoff[(abs(x.m$X) >= effect_cutoff)] <- "Significant Effect"
		}

		plot.elements <- plot.elements + geom_point(aes(X, Y, color = Cutoff), data=x.m)

		if(!missing(color)){ ## this is handled here because we want different behaviour based on if we have significance based coloring or not
			plot.elements <- plot.elements + scale_colour_manual(values = color)
		}

	} else {

		if(missing(color)){ ## this is handled here because we want different behaviour based on if we have significance based coloring or not
			color <- "black"
		}

		x.m$Cutoff <- NA_character_

		plot.elements <- plot.elements + geom_point(aes(X, Y), color = color, data=x.m)
	}


	plot.elements
}

#  Plots a PharmacoSig object into a Volcano Plot
# 
# 
# @S3method plot PharmacoSig
setMethod("plot", "PharmacoSig", plot.PharmacoSig)



