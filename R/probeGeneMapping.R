########################
## Benjamin Haibe-Kains
## All rights Reserved
## September 1, 2013
########################

#' Identify the best probes to use for each RNA gene, and reduce the data to
#' gene level
#' 
#' This function allows the summarization of RNA microarray data by the probe
#' most suitable for each gene. This function uses either the jetset package, or
#' probe variance to select the best probe to keep for each ENTREZ gene id, and
#' returns the dataset at a gene level. The different platforms are identified
#' using their GEO id.
#' 
#' @references 
#'  Li Q, Birkbak NJ, Gyorffy B, Szallasi Z and Eklund AC (2011). "Jetset:
#'  selecting the optimal microarray probe set to represent a gene." BMC
#'  Bioinformatics, 12, pp. 474.
#' 
#' @examples 
#' data(CGPsmall)
#' CGPsmall <- probeGeneMapping(CGPsmall, platform = 'GPL96', method='jetset')
#'
#' @param pSet [PharmacoSet] The PharmacoSet on which to preform the mapping platform
#' @param platform [character] identifier of the microarray platform (usually a GPL id from GEO)
#' @param method [character] either the most variant probes or jetset for Affymetrix platform
#' @return [PharmacoSet] An updated PharmacoSet object with single probe per Entrez gene id
#' @export
#' @import jetset


probeGeneMapping <- function (pSet, platform=c("MISC", "GPL8300", "GPL96", "GPL3921", "GPL97", "GPL570", "GPL1352"), method=c("variance", "jetset")){

  eset <- pSet@molecularData$rna
  
  platform <- match.arg(platform)
  method <- match.arg(method)
  
  platf.map <- rbind(c("MISC", "variance", ""),
    c("GPL8300", "jetset", "hgu95av2"),
    c("GPL96", "jetset", "hgu133a"),
    c("GPL3921", "jetset", "hgu133a"),
    c("GPL97", "jetset", "hgu133plus2"),
    c("GPL570", "jetset", "hgu133plus2"),
    c("GPL1352", "jetset", "u133x3p"))
  dimnames(platf.map) <- list(platf.map[ , 1], c("platform", "method", "parameters"))
  if (!is.element(method, platf.map[platf.map[ , "platform"], "method"])) {
    stop(sprintf("Method %s cannot be applied on platform %s\nUse the following method(s) instead: %s", method, platform, paste(x=platf.map[platf.map[ , "platform"] == platform, "method"], collapse=", ")))
  }
  params <- platf.map[which(platform == platf.map[ , "platform"]), "parameters"]
  
  
  
  switch (method,
    "jetset" = {
      
      ## keep only ENTREZID and SYMBOL in feature annotation
      # Biobase::fData(eset) <- Biobase::fData(eset)[ , c("ENTREZID", "SYMBOL"), drop=FALSE]
      # Biobase::fData(eset)[ , "ENTREZID"] <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "ENTREZID"]))
      # Biobase::fData(eset)[ , "SYMBOL"] <- stripWhiteSpace(as.character(Biobase::fData(eset)[ , "SYMBOL"]))
      #     
      js <- jetset::jscores(chip=params, probeset=rownames(Biobase::exprs(eset)))
      js <- js[rownames(Biobase::exprs(eset)), , drop=FALSE]
      ## identify the best probeset for each Entrez Gene ID
      geneid1 <- stripWhiteSpace(as.character(js[ ,"EntrezID"]))
      names(geneid1) <- rownames(js)
      geneid2 <- sort(unique(geneid1))
      names(geneid2) <- paste("geneid", geneid2, sep=".")
      gix1 <- !is.na(geneid1)
      gix2 <- !is.na(geneid2)
      geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
      ## probes corresponding to common gene ids
      gg <- names(geneid1)[is.element(geneid1, geneid.common)]
      gid <- geneid1[is.element(geneid1, geneid.common)]
      ## duplicated gene ids
      gid.dupl <- unique(gid[duplicated(gid)])
      gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
      ## unique gene ids
      gid.uniq <- gid[!is.element(gid, gid.dupl)]
      gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
      ## which are the best probe for each gene
      js <- data.frame(js, "best"=FALSE, stringsAsFactors=FALSE)
      js[gg.uniq, "best"] <- TRUE
      ## data for duplicated gene ids
      if(length(gid.dupl) > 0) {	
      	## use jetset oevrall score to select the best probesets
      	myscore <- js[gg.dupl,"overall"]
      	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
      	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
      	js[myscore[ ,"probe"], "best"] <- TRUE
      }
      ## update the esets
      probes <- rownames(Biobase::exprs(eset))[js[ , "best"]]
      names(probes) <- paste("geneid", js[js[ , "best"], "EntrezID"], sep=".")
      gid <- js[js[ , "best"], "EntrezID"]
      gsymb <- js[js[ , "best"], "symbol"]
      eset <- eset[probes, , drop=FALSE]
      rownames(Biobase::exprs(eset)) <- names(probes)
      rownames(Biobase::fData(eset)) <- names(probes)
      Biobase::fData(eset)[ , "PROBE"] <- probes
      Biobase::fData(eset)[ , "ENTREZID"] <- gid
      Biobase::fData(eset)[ , "SYMBOL"] <- gsymb
      Biobase::fData(eset)[ , "GENEID"] <- gid
    },
    "variance" = {
      ## other platform, select the most variant probe per Entrez Gene ID
      if (!all(c("GENEID", "PROBE") %in% colnames(mapping))){
        stop("Please provide columns mapping PROBE to GENEID, with those respective names")
      }
      if(!all(complete.cases(mapping))){
        warning("Only Mapped Probes were kept")
        mapping <- mapping[complete.cases(mapping),]
      }
      gid <- stripWhiteSpace(as.character(mapping[ , "GENEID"]))
      names(gid) <- as.character(mapping[ , "PROBE"])
      ugid <- sort(unique(gid))
      rr <- geneid.map(geneid1=gid, data1=t(Biobase::exprs(eset)), geneid2=ugid)
      probes <- colnames(rr$data1)
      names(probes) <- paste("geneid", rr$geneid1, sep=".")
      eset <- eset[probes, ]
      rownames(Biobase::exprs(eset)) <- names(probes)
      rownames(Biobase::fData(eset)) <- names(probes)
      Biobase::fData(eset)[ , "PROBE"] <- probes
      Biobase::fData(eset)[ , "GENEID"] <- rr$geneid1
    },
    {
      stop(sprintf("Unknown method for probe-gene mapping for platform %s", platform))
    }
  )
  pSet@molecularData$rna <- eset
  return (pSet)
  }

