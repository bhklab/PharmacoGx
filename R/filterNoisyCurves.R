filterNoisyCurves <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
    
    acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
        #for(xp in rownames(sensitivityInfo(pSet))){
        drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
        drug.responses <- drug.responses[complete.cases(drug.responses), ]
        doses.no <- nrow(drug.responses)
        
        drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
        
        delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
        
        max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
        
        if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
        (delta.sum < epsilon) &
        (max.cum.sum < (2 * epsilon)) &
        (mean(drug.responses$Viability) < mean.viablity)) {
            return (xp)
        }
    }, mc.cores=nthread)
    acceptable <- unlist(acceptable)
    noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
    return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
    xx <- as.numeric(xx)
    if(trunc)
    {
        return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
    }else{
        return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
    }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
    xx <- as.numeric(xx)
    if(trunc) {
        xx <- pmin(xx, 100)
    }
    tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
    tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
    cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
    return(max(cum.sum))
}