#### Copyright Benjamin Haibe Kains, adapted from the genefu package

`rescale` <-
function(x, na.rm=FALSE, q=0) {
	if(q == 0) {
		ma <- max(x, na.rm=na.rm)
		mi <- min(x, na.rm=na.rm)
	} else {
		ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
		mi <- quantile(x, probs=q/2, na.rm=na.rm)
	}
	xx <- (x - mi) / (ma - mi)
	return(xx)
}