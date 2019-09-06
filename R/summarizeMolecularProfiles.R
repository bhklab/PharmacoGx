#' Takes molecular data from a PharmacoSet, and summarises them
#' into one entry per drug
#'
#' Given a PharmacoSet with molecular data, this function will summarize
#' the data into one profile per cell line, using the chosed summary.stat. Note
#' that this does not really make sense with perturbation type data, and will
#' combine experiments and controls when doing the summary if run on a
#' perturbation dataset.
#'
#' @examples
#' data(GDSCsmall)
#' GDSCsmall <- summarizeMolecularProfiles(GDSCsmall,
#'                     mDataType = "rna", cell.lines=cellNames(GDSCsmall),
#'                     summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
#' GDSCsmall
#'
#' @param pSet \code{PharmacoSet} The PharmacoSet to summarize
#' @param mDataType \code{character} which one of the molecular data types
#' to use in the analysis, out of all the molecular data types available for the pset
#' for example: rna, rnaseq, snp
#' @param cell.lines \code{character} The cell lines to be summarized.
#'   If any cell.line has no data, missing values will be created
#' @param features \code{caracter} A vector of the feature names to include in the summary
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell.lines? Choices are "mean", "median", "first", or "last"
#'   In case molecular data type is mutation or fusion "and" and "or" choices are available 
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param summarize A flag which when set to FALSE (defaults to TRUE) disables summarizing and
#'   returns the data unchanged as a ExpressionSet
#' @param verbose \code{boolean} should messages be printed
#' @param binarize.threshold \code{numeric} A value on which the molecular data is binarized. 
#'   If NA, no binarization is done.
#' @param binarize.direction \code{character} One of "less" or "greater", the direction of binarization on 
#'   binarize.threshold, if it is not NA. 
#' @return \code{matrix} An updated PharmacoSet with the molecular data summarized
#'   per cell line.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<- colData colData<- assays assays<-
#' @importFrom Biobase AnnotatedDataFrame
#' @export

##TODO:: Add features parameter

summarizeMolecularProfiles <- function(pSet, 
                                       mDataType, 
                                       cell.lines, 
                                       features, 
                                       summary.stat = c("mean", "median", "first", "last", "and", "or"), 
                                       fill.missing = TRUE, 
                                       summarize = TRUE, 
                                       verbose = TRUE,
                                       binarize.threshold = NA, 
                                       binarize.direction = c("less", "greater")
                                       ) {
  
  
  ### Placed here to make sure the pSet argument gets checked first by R. 
  mDataTypes <- names(pSet@molecularProfiles)
  if (!(mDataType %in% mDataTypes)) {
    stop (sprintf("Invalid mDataType, choose among: %s", paste(names(pSet@molecularProfiles), collapse=", ")))
  }
  
  if(summarize==FALSE){
    return(pSet@molecularProfiles[[mDataType]])
  }
  
  if (missing(features)) {
    features <- rownames(featureInfo(pSet, mDataType))
  } else {
    fix <- is.element(features, rownames(featureInfo(pSet, mDataType)))
    if (verbose && !all(fix)) {
      warning (sprintf("Only %i/%i features can be found", sum(fix), length(features)))
    }
    features <- features[fix]
  }
  
  summary.stat <- match.arg(summary.stat)
  binarize.direction <- match.arg(binarize.direction)
  
  if((!S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation %in% c("mutation","fusion")) & (!summary.stat %in% c("mean", "median", "first", "last"))) {
    stop ("Invalid summary.stat, choose among: mean, median, first, last" )
  }
  if((S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation %in% c("mutation","fusion")) & (!summary.stat %in% c("and", "or"))) {
    stop ("Invalid summary.stat, choose among: and, or" )
  }
  
  if (missing(cell.lines)) {
    cell.lines <- cellNames(pSet)
  }
  
  dd <- molecularProfiles(pSet, mDataType)
  pp <- phenoInfo(pSet, mDataType)
  
  if(S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation == "mutation") {
    tt <- dd
    tt[which(!is.na(dd) & dd =="wt")] <- FALSE
    tt[which(!is.na(dd) & dd !="wt")] <- TRUE
    tt <- apply(tt, 2, as.logical)
    dimnames(tt) <- dimnames(dd)
    dd <- tt
  }
  if(S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation == "fusion") {
    tt <- dd
    tt[which(!is.na(dd) & dd =="0")] <- FALSE
    tt[which(!is.na(dd) & dd !="0")] <- TRUE
    tt <- apply(tt, 2, as.logical)
    dimnames(tt) <- dimnames(dd)
    dd <- tt
  }
  if(S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation %in% c("cnv", "rna", "rnaseq", "isoform") 
     && !is.na(binarize.threshold)) {
    tt <- dd
    switch(binarize.direction, "less" = {
          tt[which(!is.na(dd) & dd < binarize.threshold)] <- TRUE
          tt[which(!is.na(dd) & dd >= binarize.threshold)] <- FALSE
    }, "greater" = {
          tt[which(!is.na(dd) & dd > binarize.threshold)] <- TRUE
          tt[which(!is.na(dd) & dd <= binarize.threshold)] <- FALSE
    })
    tt <- apply(tt, 2, as.logical)
    dimnames(tt) <- dimnames(dd)
    dd <- tt
  }
  if (any(colnames(dd) != rownames(pp))) {
    warning ("Samples in phenodata and expression matrices must be ordered the same way")
    dd <- dd[ , rownames(pp), drop=FALSE]
  }
  if (!fill.missing) {
    cell.lines <- intersect(cell.lines, unique(pp[!is.na(pp[ , "cellid"]), "cellid"]))
  }
  if (length(cell.lines) == 0) {
    stop ("No cell lines in common")
  }
    
  ## select profiles with no replicates
  duplix <- unique(pp[!is.na(pp[ , "cellid"]) & duplicated(pp[ , "cellid"]), "cellid"])
  ucell <- setdiff(cell.lines, duplix)
  
  ## keep the non ambiguous cases
  dd2 <- dd[ , match(ucell, pp[ , "cellid"]), drop=FALSE]
  pp2 <- pp[match(ucell, pp[ , "cellid"]), , drop=FALSE]
  if (length(duplix) > 0) {
    if (verbose) {
      message(sprintf("Summarizing %s molecular data for:\t%s", mDataType, pSet@annotation$name))
      total <- length(duplix)
      # create progress bar 
      pb <- utils::txtProgressBar(min=0, max=total, style=3)
      i <- 1
    }
    ## replace factors by characters to allow for merging duplicated experiments
    pp2 <- apply(pp2, 2, function (x) {
      if (is.factor(x)) {
        return (as.character(x))
      } else { 
        return (x) 
      } 
    })
    ## there are some replicates to collapse
    for (x in duplix) {
      myx <- which(!is.na(pp[ , "cellid"]) & is.element(pp[ , "cellid"], x))
      switch(summary.stat,
        "mean" = {
          ddt <- apply(dd[ , myx, drop=FALSE], 1, mean)
        },
        "median"={
          ddt <- apply(dd[ , myx, drop=FALSE], 1, median)
        }, 
        "first"={
          ddt <- dd[ , myx[1], drop=FALSE]
        },
        "last" = {
          ddt <- dd[ , myx[length(myx)], drop=FALSE]
        },
        "and" = {
          ddt <- apply(dd[ , myx, drop=FALSE], 1, function(x) do.call(`&`, as.list(x)))
        },
        "or" = {
          ddt <- apply(dd[ , myx, drop=FALSE], 1, function(x) do.call(`|`, as.list(x)))
        }
      )
      ppt <- apply(pp[myx, , drop=FALSE], 2, function (x) {
        x <- paste(unique(as.character(x[!is.na(x)])), collapse="///")
        return (x)
      })
      ppt[!is.na(ppt) & ppt == ""] <- NA
      dd2 <- cbind(dd2, ddt)
      pp2 <- rbind(pp2, ppt)
      if (verbose){
        utils::setTxtProgressBar(pb, i)
        i <- i + 1
      }
    }
    if (verbose) {
      close(pb)
    }
  }
  colnames(dd2) <- rownames(pp2) <- c(ucell, duplix)
  
  ## reorder cell lines
  dd2 <- dd2[ , cell.lines, drop=FALSE]
  pp2 <- pp2[cell.lines, , drop=FALSE]
  pp2[ , "cellid"] <- cell.lines
  res <- pSet@molecularProfiles[[mDataType]]
  if(S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])$annotation %in% c("mutation", "fusion")) {
    tt <- dd2
    tt[which(!is.na(dd2) & dd2)] <- "1"
    tt[which(!is.na(dd2) & !dd2)] <- "0"
    dd2 <- tt
  }
  res <- SummarizedExperiment::SummarizedExperiment(dd2)
  #Biobase::exprs(res) <- dd2
  pp2 <- S4Vectors::DataFrame(pp2, row.names=rownames(pp2))
  pp2$tissueid <- cellInfo(pSet)[pp2$cellid, "tissueid"]
  SummarizedExperiment::colData(res) <- pp2
  SummarizedExperiment::rowData(res) <- featureInfo(pSet, mDataType)
  #Biobase::exprs(res) <- Biobase::exprs(res)[features,]
  #Biobase::fData(res) <- Biobase::fData(res)[features,]
  res <- res[features,]
  S4Vectors::metadata(res) <- S4Vectors::metadata(pSet@molecularProfiles[[mDataType]])
  ####NOTE:: Removed se.exprs from assays in conversion to SE
  #if(!is.null(SummarizedExperiment::assays(res)$exprs)) SummarizedExperiment::assays(res)$se.exprs <- NULL
  return(res)
}
