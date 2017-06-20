#' Intersects objects of the PharmacoSet class, subsetting them to the common
#' drugs and/or cell lines as selected by the user.
#' 
#' Given a list of PharmacoSets, the function will find the common drugs,
#' and/or cell lines, and return PharmacoSets that contain data only pertaining 
#' to the common drugs, and/or cell lines. The mapping between dataset 
#' drug and cell names is done using annotations found in the 
#' PharmacoSet object's internal curation slot
#' 
#' @examples 
#' data(GDSCsmall)
#' data(CCLEsmall)
#' common <- intersectPSet(list('GDSC'=GDSCsmall,
#'   'CCLE'=CCLEsmall), intersectOn = c("drugs", "cell.lines", "concentrations"))
#' common$CGP
#' common$CCLE
#' 
#' 
#' @param pSets [list] a list of PharmacoSet objects, of which the function
#'   should find the intersection
#' @param intersectOn [character] which identifiers to intersect on, 
#'   drugs or cell lines
#' @param drugs a vector of common drugs between pSets.
#' In case user is intersted on getting intersection on certain drugs,
#' they can provide their list of drugs.
#' @param cells a vector of common cell lines between pSets.
#' In case user is intersted on getting intersection on certain cell lines,
#' they can provide their list of cell lines
#' @param strictIntersect [boolean] Should the intersection keep only the drugs 
#'   and cell lines that have been tested on together?
#' @param verbose [boolean] Should the function announce its key steps?
#' @param nthread [numeric] The number of cores to use to run intersection on concentrations
#' @return [list] a list of pSets, contatining only the intersection
#' @export
#' 
intersectPSet <- function (pSets, intersectOn=c("drugs", "cell.lines", "concentrations"), cells, drugs, strictIntersect=FALSE, verbose=TRUE, nthread=1) {
  
  if (verbose) {
    message("Intersecting large PSets may take a long time ...")
  }
  ## TODO: Fix the strict intersection!!!!!!
  if (length(pSets) == 1) {
    return(pSets) 
  }
  if (length(pSets) > 1) {
    if(is.null(names(pSets)) ){
      
      names(pSets) <- sapply(pSets, pSetName)
      
    }
    if ("drugs" %in% intersectOn){
      common.drugs <- intersectList(lapply(pSets, function(x){return(drugNames(x))}))
      if(!missing(drugs)) {
        common.drugs <- intersect(common.drugs, drugs)
      }
      if (length(common.drugs) == 0) {
        stop("No drugs is in common between pSets!")
      }
    }
    if ("cell.lines" %in% intersectOn){
      common.cells <- intersectList(lapply(pSets, function(x){return(cellNames(x))}))
      if(!missing(cells)) {
        common.cells <- intersect(common.cells, cells)
      }
      if (length(common.cells) == 0) {
        stop("No cell lines is in common between pSets!")
      }
    } 
    if (("drugs" %in% intersectOn) & ("cell.lines" %in% intersectOn)) {
      common.exps <- intersectList(lapply(pSets, function (x){
        if ("cellid" %in% colnames(x@sensitivity$info) & "drugid" %in% colnames(x@sensitivity$info)) {
          paste(x@sensitivity$info$cellid, x@sensitivity$info$drugid, sep = "_")
        } else { NULL }
      }))
      # expMatch <- data.frame(lapply(pSets, 
      #   function (x, common.exps){
      #     if ("cellid" %in% colnames(x@sensitivity$info) & "drugid" %in% colnames(x@sensitivity$info)){
      
      #       myx <- match(paste(x@sensitivity$info$cellid, x@sensitivity$info$drugid, sep = "_") ,common.exps)
      
      #       res <- rownames(x@sensitivity$info)[!is.na(myx)]
      
      #       names(res) <- common.exps[na.omit(myx)]
      
      #       res <- res[common.exps]
      
      #       return(res)
      
      #     } else { NULL }
      #   }, common.exps=common.exps))
      expMatch <- lapply(pSets, 
                         function (x, common.exps){
                           if ("cellid" %in% colnames(x@sensitivity$info) & "drugid" %in% colnames(x@sensitivity$info)){
                             
                             myx <- match(paste(x@sensitivity$info$cellid, x@sensitivity$info$drugid, sep = "_") ,common.exps)
                             
                             res <- rownames(x@sensitivity$info)[!is.na(myx)]
                             
                             names(res) <- common.exps[na.omit(myx)]
                             
                             res <- res[common.exps]
                             
                             return(res)
                             
                           } else { NULL }
                         }, common.exps=common.exps)
      # }, common.exps=common.exps)
      
      if(strictIntersect){
        if(length(unique(sapply(expMatch, length)))>1){
          stop("Strict Intersecting works only when each PSet has 1 replicate per cell-drug pair. Use collapseSensitvityReplicates to reduce the sensitivity data as required")
        }
        expMatch <- data.frame(expMatch,  stringsAsFactors=FALSE)
        # expMatch2 <- as.matrix(expMatch2)
        rownames(expMatch) <- common.exps
        colnames(expMatch) <- names(pSets)
        
      } else {
        
        expMatch <- lapply(expMatch, function(x){names(x) <- x; return(x)})
      }
    }
    if (("drugs" %in% intersectOn) & ("cell.lines" %in% intersectOn) & ("concentrations" %in% intersectOn)) {
      
      if(length(unique(sapply(expMatch, length)))>1){
        stop("Intersecting on concentrations works only when each PSet has 1 replicate per cell-drug pair. Use collapseSensitvityReplicates to reduce the sensitivity data as required")
      }
      
      expMatch <- data.frame(expMatch,  stringsAsFactors=FALSE)
      # expMatch2 <- as.matrix(expMatch2)
      rownames(expMatch) <- common.exps
      colnames(expMatch) <- names(pSets)
      
      pSets <- .calculateSensitivitiesStar(pSets, exps=expMatch, cap=100, nthread=nthread)
    }
    if ("cell.lines" %in% intersectOn){
      molecular.types  <- NULL
      for (pSet in pSets) {
        for (eSet in pSet@molecularProfiles) {
          molecular.types <- union(molecular.types, ifelse(length(grep("rna", Biobase::annotation(eSet)) > 0), "rna", Biobase::annotation(eSet)))
        }
      }
      common.molecular.cells <- list()
      for (molecular.type in molecular.types) {
        if(strictIntersect){
          common.molecular.cells[[molecular.type]] <- intersectList(lapply(pSets, 
                                                                           function (pSet) { 
                                                                             eSets <- names(unlist(sapply(pSet@molecularProfiles, function(eSet){grep(molecular.type, Biobase::annotation(eSet))})))
                                                                             if(length(eSets) > 0){
                                                                               return(intersectList(sapply(eSets, 
                                                                                                           function(eSet) {
                                                                                                             if (length(grep(molecular.type, Biobase::annotation(pSet@molecularProfiles[[eSet]]))) > 0) {
                                                                                                               intersect(Biobase::pData(pSet@molecularProfiles[[eSet]])$cellid, common.cells)
                                                                                                             }
                                                                                                           })))
                                                                             }
                                                                           }))
        }else{
          common.molecular.cells[[molecular.type]] <- intersectList(lapply(pSets, 
                                                                           function (pSet) { 
                                                                             eSets <- names(unlist(sapply(pSet@molecularProfiles, function(eSet){grep(molecular.type, Biobase::annotation(eSet))})))
                                                                             return(unionList(sapply(eSets, 
                                                                                                     function(eSet) {
                                                                                                       if (length(grep(molecular.type, Biobase::annotation(pSet@molecularProfiles[[eSet]]))) > 0) {
                                                                                                         intersect(Biobase::pData(pSet@molecularProfiles[[eSet]])$cellid, common.cells)
                                                                                                       }
                                                                                                     })))
                                                                           }))
        }
      }
      
      
    }
    
    
    
    
    for(i in 1:length(pSets)){
      if(("drugs" %in% intersectOn) & ("cell.lines" %in% intersectOn)){
        if(strictIntersect){
          pSets[[i]] <- subsetTo(pSet=pSets[[i]], drugs=common.drugs, cells=common.cells, exps=expMatch, molecular.data.cells=common.molecular.cells)
          
        } else {
          pSets[[i]] <- subsetTo(pSet=pSets[[i]], drugs=common.drugs, cells=common.cells, molecular.data.cells=common.molecular.cells)
        } 
      } else if(("cell.lines" %in% intersectOn)) {
        pSets[[i]] <- subsetTo(pSet=pSets[[i]], cells=common.cells, molecular.data.cells=common.molecular.cells)
        
      } else if(("drugs" %in% intersectOn)) {
        pSets[[i]] <- subsetTo(pSet=pSets[[i]], drugs=common.drugs)
        
      }
    }
    return(pSets)
  }
}


