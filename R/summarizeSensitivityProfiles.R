#' Takes the sensitivity data from a PharmacoSet, and summarises them into a
#' drug vs cell line table
#' 
#' This function creates a table with cell lines as rows and drugs as columns,
#' summarising the drug senstitivity data of a PharmacoSet into drug-cell line
#' pairs
#' 
#' @examples 
#' data(GDSCsmall)
#' GDSCauc <- summarizeSensitivityProfiles(GDSCsmall, sensitivity.measure='auc_published')
#'
#' @param pSet [PharmacoSet] The PharmacoSet from which to extract the data
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? The current
#'   choices are 'ic50_published', 'auc_published', 'ic50_recomputed',
#'   'auc_recomputed', 'auc_recomputed_star'.
#' @param cell.lines \code{character} The cell lines to be summarized. 
#'    If any cell lines has no data, it will be filled with
#'   missing values
#' @param drugs \code{character} The drugs to be summarized.
#'   If any drugs has no data, it will be filled with
#'   missing values
#' @param summary.stat \code{character} which summary method to use if there are repeated
#'   cell line-drug experiments? Choices are "mean", "median", "first", or "last"
#' @param fill.missing \code{boolean} should the missing cell lines not in the
#'   molecular data object be filled in with missing values?
#' @param verbose Should the function print progress messages?
#' @return [matrix] A matrix with cell lines going down the rows, drugs across
#'   the columns, with the selected sensitivity statistic for each pair.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats median
#' @export


summarizeSensitivityProfiles <- function(pSet, sensitivity.measure=c("ic50_published", "auc_published", "ic50_recomputed", "auc_recomputed", "auc_recomputed_star", "Synergy_score", "amax_published", "amax_recomputed"), cell.lines, drugs, summary.stat=c("mean", "median", "first", "last"), fill.missing=TRUE, verbose=TRUE){
    
    summary.stat <- match.arg(summary.stat)
  sensitivity.measure <- match.arg(sensitivity.measure)
  if (!(sensitivity.measure %in% colnames(sensitivityProfiles(pSet)))) {
    stop (sprintf("Invalid sensitivity measure for %s, choose among: %s", pSet@annotation$name, paste(colnames(sensitivityProfiles(pSet)), collapse=", ")))
  }
  if (missing(cell.lines)) {
    cell.lines <- cellNames(pSet)
  }
  if (missing(drugs)) {
    if(sensitivity.measure != "Synergy_score")
    {
      drugs <- drugNames(pSet)
    }else{
      drugs <- sensitivityInfo(pSet)[grep("///", sensitivityInfo(pSet)$drugid), "drugid"]
    }
  }
  
    pp <- sensitivityInfo(pSet)
  pp <- pp[which(pp$cellid %in% cell.lines & pp$drugid %in% drugs),]
    dd <- sensitivityProfiles(pSet)[rownames(pp),]
    
    if (!fill.missing) {
      cell.lines <- intersect(cell.lines, unique(pp[!is.na(pp[ , "cellid"]), "cellid"]))
    }
    if (!fill.missing) {
      drugs <- intersect(drugs, unique(pp[!is.na(pp[ , "drugid"]), "drugid"]))
    }
    
    
    ## select profiles with no replicates
  # xps <- apply(pp[ , c("drugid", "cellid")], 1, function (x) {
  #   if(any(is.na(x))) {
  #     x <- NA
  #   } else {
  #     x <- paste(x, collapse="///")
  #   }
  #   return (x)
  # })
  # names(xps) <- rownames(pp)
  
  xps <- apply(pp[ , c("drugid", "cellid")], 1, paste, collapse="////")
  xps[!complete.cases(pp[ , c("drugid", "cellid")])] <- NA
  duplix <- unique(xps[!is.na(xps) & duplicated(xps)])
  uniqix <- setdiff(xps[!is.na(xps)], duplix)
  iix <- t(sapply(strsplit(uniqix, "////"), function (x) { return (x) }))
  iix <- cbind(match(iix[ , 1], drugs), match(iix[ , 2], cell.lines))
  iix2 <- match(uniqix, xps)
  ## keep the non ambiguous cases
  dd2 <- matrix(NA, nrow=length(drugs), ncol=length(cell.lines), dimnames=list(drugs, cell.lines))
  for (ii in 1:nrow(iix)) {
    dd2[iix[ii, 1], iix[ii, 2]] <- dd[iix2[ii], sensitivity.measure]
  }
  # pp2 <- pp[match(uniqix, xps), , drop=FALSE]
  # rownames(pp2) <- uniqix
  
  if (length(duplix) > 0) {
    if (verbose) {
      message(sprintf("Summarizing %s sensitivity data for:\t%s", sensitivity.measure, pSet@annotation$name))
      total <- length(duplix)
      # create progress bar 
      pb <- utils::txtProgressBar(min=0, max=total, style=3)
      i <- 1
    }
    ## there are some replicates to collapse
    for (x in duplix) {
      myx <- which(!is.na(xps) & xps == x)
      iix <- unlist(strsplit(duplix, "////"))
            switch(summary.stat, 
        "mean" = {
                  dd2[iix[1], iix[2]] <- mean(dd[myx, sensitivity.measure])
              },
              "median" = {
                  dd2[iix[1], iix[2]] <- median(dd[myx, sensitivity.measure])
              }, 
              "first" = {
                  dd2[iix[1], iix[2]] <- dd[myx[1], sensitivity.measure]
              },
              "last" = {
                  dd2[iix[1], iix[2]] <- dd[myx[length(myx)], sensitivity.measure]
              }
      )
      # ppt <- apply(pp[myx, , drop=FALSE], 2, function (x) {
      #   x <- paste(unique(x), collapse="////")
      #   return (x)
      # })
      # pp2 <- rbind(pp2, ppt)
      if (verbose){
        utils::setTxtProgressBar(pb, i)
        i <- i + 1
      }
    }
    if (verbose) {
      close(pb)
    }
  }
  res <- dd2
  ## TODO: return the collapsed sensitivty Info as well?
  if(sensitivity.measure != "Synergy_score") {
    return(res)
  }else{
    if(fill.missing){
      dd <- drugNames(pSet)
    } else {
      dd <- unique(unlist(strsplit(drugs, split="///")))
    }
    tt <- array(NA, dim=c(length(dd), length(dd), length(cell.lines)), dimnames=list(dd, dd, cell.lines))
    for(drug in dd) {
      tt[drug, dd, cell.lines] <- res[match(sprintf("%s///%s", drug, dd), rownames(res)), cell.lines]
    }
    return(tt)
  }
}