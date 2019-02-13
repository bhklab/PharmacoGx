######
# Checks known biomakers given a mCI style signature file.
#####


#' Check a mCI style signature file for the results of known biomarkers
#' 
#' @param res mCI style drug Sensitivity Signature result
#' 
#' @export


checkKnownBiomarkers <- function(res, type="expression"){

	to.check <- subset(PharmacoGx::known.biomarkers, Type == type)

  if(type == "expression"){
    to.check <- to.check[!is.na(to.check$Ensembl), ]
  } 

  pan.cancer <- res[[1]]
  tissue.specific <- res[[3]]



  my.genes <- as.character(to.check$Ensembl[to.check$Ensembl %in% dimnames(pan.cancer)[[3]]])
  my.drugs <- as.character(to.check$Drug[to.check$Drug %in% dimnames(pan.cancer)[[1]]])

  known.res <- as.data.frame(t(mapply(function(x,y){
    return(pan.cancer[x,,y])}, my.drugs, my.genes)))

  known.res <- data.frame("Ensembl" = my.genes, "Drug" = my.drugs, 
                          "Symbol" = to.check$Gene[to.check$Ensembl %in% dimnames(pan.cancer)[[3]]],
                          known.res)

  known.res.tissue <- as.data.frame(t(mapply(function(x,y){
    return(tissue.specific[x,"CI",y,])}, my.drugs, my.genes)))
  known.res.tissue.p <- as.data.frame(t(mapply(function(x,y){
    return(tissue.specific[x,"fdr",y,])}, my.drugs, my.genes)))
  return(list(known.res, known.res.tissue, known.res.tissue.p))

}

require(ggplot2)
investigateKnownBiomarkers <- function(res, pSet, type="expression", mDataType="rnaseq") {

  known.list <- checkKnownBiomarkers(res)
  known.res <- known.list[[1]]
  known.res.tissue <- known.list[[2]]
  known.res.tissue.p <- known.list[[3]] 

  molecularProfiles <- Biobase::exprs(summarizeMolecularProfiles(pSet, mDataType))
  drugSensitivity  <- summarizeSensitivityProfiles(pSet)

  tissue_ids <- cellInfo(pSet)$tissueid

  pdf(paste0("known_biomarker_report_", pSetName(pSet), "_mDataType_", mDataType, ".pdf"))
  for(i in 1:NROW(known.res)){

    toPlot <- data.frame(molecularProfiles[known.res$Ensembl[i],],
                         drugSensitivity[known.res$Drug[i],],
                         tissue_ids)
    colnames(toPlot) <- make.names(c(as.character(known.res$Symbol[i]), as.character(known.res$Drug[i]), "tissue"))


    anotCI <- known.res[i,"CI", drop=FALSE]
    anotFDR <- known.res[i,"fdr", drop=FALSE]
    p <- ggplot(toPlot, aes_string(x = colnames(toPlot)[1], y = colnames(toPlot)[2]), colour = "tissue") + scale_x_log10() +geom_point(aes(colour = tissue)) + theme_classic() + ggtitle("Pan-Cancer") + theme(legend.position="none") + geom_label(data = anotCI, aes(label=CI), x=max(log10(toPlot[[1]]), na.rm=TRUE), y =max(toPlot[[2]], na.rm=TRUE), hjust=1, vjust=1) + geom_label(data = anotFDR, aes(label=fdr), x=max(log10(toPlot[[1]]), na.rm=TRUE), y =max(toPlot[[2]], na.rm=TRUE), hjust=1, vjust=2)
    print(p)
    # grid::grid.newpage()
  }
  for(i in 1:NROW(known.res)){

    toPlot <- data.frame(molecularProfiles[as.character(known.res$Ensembl[i]),],
                         drugSensitivity[as.character(known.res$Drug[i]),],
                         tissue_ids)
    colnames(toPlot) <- make.names(c(as.character(known.res$Symbol[i]), as.character(known.res$Drug[i]), "tissue"))

    anotCI <- data.frame("tissue" = na.omit(unique(tissue_ids)), "CI" = round(as.numeric(known.res.tissue[i,na.omit(unique(tissue_ids))]),2))
    anotFDR <- data.frame("tissue" = na.omit(unique(tissue_ids)), "p" = round(as.numeric(known.res.tissue.p[i,na.omit(unique(tissue_ids))]),2))

    p <- ggplot(toPlot, aes_string(x = colnames(toPlot)[1], y = colnames(toPlot)[2])) + scale_x_log10() + geom_point() + theme_classic() + ggtitle("Pan-Cancer") + theme(legend.position="none") + facet_wrap(~ tissue, ncol = 5) + geom_label(data = anotCI, aes(label=CI), x=max(log10(toPlot[[1]]), na.rm=TRUE), y =max(toPlot[[2]], na.rm=TRUE), hjust=1, vjust=1) + geom_label(data = anotFDR, aes(label=p), x=max(log10(toPlot[[1]]), na.rm=TRUE), y =max(toPlot[[2]], na.rm=TRUE), hjust=1, vjust=2)
    
    print(p)
    # grid::grid.newpage()
  }
  dev.off()

}