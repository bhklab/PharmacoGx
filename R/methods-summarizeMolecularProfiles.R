#' Takes molecular data from a PharmacoSet, and summarises them
#' into one entry per drug
#'
#' Given a PharmacoSet with molecular data, this function will summarize
#' the data into one profile per cell line, using the chosen summary.stat. Note
#' that this does not really make sense with perturbation type data, and will
#' combine experiments and controls when doing the summary if run on a
#' perturbation dataset.
#'
#' @examples
#' data(GDSCsmall)
#' GDSCsmall <- summarizeMolecularProfiles(GDSCsmall, mDataType = "rna", cell.lines=sampleNames(GDSCsmall), summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)
#' GDSCsmall
#'
#' @param object \code{PharmacoSet} The PharmacoSet to summarize
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
#' @param removeTreated \code{logical} If treated/perturbation experiments are present, should they
#'   be removed? Defaults to yes.
#'
#' @return \code{matrix} An updated PharmacoSet with the molecular data summarized
#'   per cell line.
#'
#' @importMethodsFrom CoreGx summarizeMolecularProfiles
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<- colData colData<- assays assays<- assayNames assayNames<-
#' @importFrom Biobase AnnotatedDataFrame
#' @keywords internal
#' @export
setMethod('summarizeMolecularProfiles', signature(object='PharmacoSet'),
    function(object, mDataType, cell.lines, features,
        summary.stat = c("mean", "median", "first", "last", "and", "or"),
        fill.missing = TRUE, summarize = TRUE, verbose = TRUE,
        binarize.threshold = NA, binarize.direction = c("less", "greater"),
        removeTreated=TRUE)
{

    mDataTypes <- mDataNames(object)
    if (!(mDataType %in% mDataTypes)) {
      stop (sprintf("Invalid mDataType, choose among: %s", paste(names(molecularProfilesSlot(object)), collapse=", ")))
    }

    if(summarize==FALSE){
      return(molecularProfilesSlot(object)[[mDataType]])
    }

    if (missing(features)) {
      features <- rownames(featureInfo(object, mDataType))
    } else {
      fix <- is.element(features, rownames(featureInfo(object, mDataType)))
      if (verbose && !all(fix)) {
        warning (sprintf("Only %i/%i features can be found", sum(fix), length(features)))
      }
      features <- features[fix]
    }

    summary.stat <- match.arg(summary.stat)
    binarize.direction <- match.arg(binarize.direction)

    if((!S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation %in% c("mutation","fusion")) & (!summary.stat %in% c("mean", "median", "first", "last"))) {
      stop ("Invalid summary.stat, choose among: mean, median, first, last" )
    }
    if((S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation %in% c("mutation","fusion")) & (!summary.stat %in% c("and", "or"))) {
      stop ("Invalid summary.stat, choose among: and, or" )
    }

    if (missing(cell.lines)) {
      cell.lines <- sampleNames(object)
    }

    if(datasetType(object) %in% c("perturbation", "both") && removeTreated){
      if(!"xptype" %in% colnames(phenoInfo(object, mDataType))) {
        warning("The passed in molecular data had no column: xptype.
                 \rEither the mDataType does not include perturbations, or the PSet is malformed.
                 \rAssuming the former and continuing.")
      } else {
        keepCols <- phenoInfo(object, mDataType)$xptype %in% c("control", "untreated")
        molecularProfilesSlot(object)[[mDataType]] <- molecularProfilesSlot(object)[[mDataType]][,keepCols]
      }
    }


    ##TODO:: have less confusing variable names
    dd <- molecularProfiles(object, mDataType)
    pp <- phenoInfo(object, mDataType)

    if(S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation == "mutation") {
      tt <- dd
      tt[which(!is.na(dd) & dd =="wt")] <- FALSE
      tt[which(!is.na(dd) & dd !="wt")] <- TRUE
      tt <- apply(tt, 2, as.logical)
      dimnames(tt) <- dimnames(dd)
      dd <- tt
    }
    if(S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation == "fusion") {
      tt <- dd
      tt[which(!is.na(dd) & dd =="0")] <- FALSE
      tt[which(!is.na(dd) & dd !="0")] <- TRUE
      tt <- apply(tt, 2, as.logical)
      dimnames(tt) <- dimnames(dd)
      dd <- tt
    }
    if(S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation %in% c("cnv", "rna", "rnaseq", "isoform")
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
      cell.lines <- intersect(cell.lines, unique(pp[!is.na(pp[ , "sampleid"]), "sampleid"]))
    }
    if (length(cell.lines) == 0) {
      stop ("No cell lines in common")
    }

    ## select profiles with no replicates
    duplix <- unique(pp[!is.na(pp[ , "sampleid"]) & duplicated(pp[ , "sampleid"]), "sampleid"])
    ucell <- setdiff(cell.lines, duplix)

    ## keep the non ambiguous cases
    dd2 <- dd[ , match(ucell, pp[ , "sampleid"]), drop=FALSE]
    pp2 <- pp[match(ucell, pp[ , "sampleid"]), , drop=FALSE]
    if (length(duplix) > 0) {
      if (verbose) {
        message(sprintf("Summarizing %s molecular data for:\t%s", mDataType, annotation(object)$name))
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
        myx <- which(!is.na(pp[ , "sampleid"]) & is.element(pp[ , "sampleid"], x))
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
    pp2[ , "sampleid"] <- cell.lines
    res <- molecularProfilesSlot(object)[[mDataType]]
    if(S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])$annotation %in% c("mutation", "fusion")) {
      tt <- dd2
      tt[which(!is.na(dd2) & dd2)] <- "1"
      tt[which(!is.na(dd2) & !dd2)] <- "0"
      dd2 <- tt
    }
    res <- SummarizedExperiment::SummarizedExperiment(dd2)
    pp2 <- S4Vectors::DataFrame(pp2, row.names=rownames(pp2))
    pp2$tissueid <- sampleInfo(object)[pp2$sampleid, "tissueid"]
    SummarizedExperiment::colData(res) <- pp2
    SummarizedExperiment::rowData(res) <- featureInfo(object, mDataType)
    ##TODO:: Generalize this to multiple assay SummarizedExperiments!
    # if(!is.null(SummarizedExperiment::assay(res, 1))) {
    #   SummarizedExperiment::assay(res, 2) <- matrix(rep(NA,
    #                                                     length(assay(res, 1))
    #                                                     ),
    #                                                     nrow=nrow(assay(res, 1)),
    #                                                     ncol=ncol(assay(res, 1)),
    #                                                 dimnames=dimnames(assay(res, 1))
    #                                                 )
    # }
    assayNames(res) <- assayNames(molecularProfilesSlot(object)[[mDataType]])[[1]]
    res <- res[features,]
    S4Vectors::metadata(res) <- S4Vectors::metadata(molecularProfilesSlot(object)[[mDataType]])
    return(res)
})
