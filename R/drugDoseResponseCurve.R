#' Plot drug response curve of a given drug and a given cell for a list of pSets (objects of the PharmacoSet class).
#' 
#' Given a list of PharmacoSets, the function will plot the drug_response curve,
#' for a given drug/cell pair. The y axis of the plot is the viability percentage
#' and x axis is the log transformed concentrations. If more than one pSet is 
#' provided, a light gray area would show the common concentration range between pSets. 
#' User can ask for type of sensitivity measurment to be shown in the plot legend.
#' The user can also provide a list of their own concentrations and viability values, 
#' as in the examples below, and it will be treated as experiments equivalent to values coming 
#' from a pset. The names of the concentration list determine the legend labels. 
#' 
#' @examples
#' if (interactive()) {
#' drugDoseResponseCurve(concentrations=list("Experiment 1"=c(.008, .04, .2, 1)),
#'  viabilities=list(c(100,50,30,1)), plot.type="Both")
#' }
#' 
#' @param drug [string] A drug name for which the drug response curve should be 
#' plotted. If the plot is desirable for more than one pharmaco set, A unique drug id
#' should be provided.
#' @param cellline [string] A cell line name for which the drug response curve should be 
#' plotted. If the plot is desirable for more than one pharmaco set, A unique cell id
#' should be provided.
#' @param pSets [list] a list of PharmacoSet objects, for which the function
#' should plot the curves.
#' @param concentrations,viabilities [list] A list of concentrations and viabilities to plot, the function assumes that
#' concentrations[[i]] is plotted against viabilities[[i]]. The names of the concentration list are used to create the legend labels
#' @param conc_as_log [logical], if true, assumes that log10-concentration data has been given rather than concentration data,
#' and that log10(ICn) should be returned instead of ICn. Applies only to the concentrations parameter.
#' @param viability_as_pct [logical], if false, assumes that viability is given as a decimal rather
#' than a percentage, and that E_inf passed in as decimal. Applies only to the viabilities parameter.
#' @param legends.label [vector] A vector of sensitivity measurment types which could 
#' be any combination of  ic50_published, auc_published, auc_recomputed and auc_recomputed_star.
#' A legend will be displayed on the top right of the plot which each line of the legend is 
#' the values of requested sensitivity measerments for one of the requested pSets.
#' If this parameter is missed no legend would be provided for the plot.
#' @param ylim [vector] A vector of two numerical values to be used as ylim of the plot.
#' If this parameter would be missed c(0,100) would be used as the ylim of the plot.
#' @param xlim [vector] A vector of two numerical values to be used as xlim of the plot.
#' If this parameter would be missed the minimum and maximum comncentrations between all 
#' the pSets would be used as plot xlim.
#' @param mycol [vector] A vector with the same lenght of the pSets parameter which 
#' will determine the color of the curve for the pharmaco sets. If this parameter is 
#' missed default colors from Rcolorbrewer package will be used as curves color. 
#' @param plot.type [character] Plot type which can be the actual one ("Actual") or 
#' the one fitted by logl logistic regression ("Fitted") or both of them ("Both").
#' If this parameter is missed by default actual curve is plotted. 
#' @param summarize.replicates [character] If this parameter is set to true replicates 
#' are summarized and replicates are plotted individually otherwise
#' @param title [character] The title of the graph. If no title is provided, then it defaults to
#' 'Drug':'Cell Line'.
#' @param lwd [numeric] The line width to plot with
#' @param cex [numeric] The cex parameter passed to plot
#' @param cex.main [numeric] The cex.main parameter passed to plot, controls the size of the titles
#' @param legend.loc And argument passable to xy.coords for the position to place the legend. 
#' @param trunc [bool] Should the viability values be truncated to lie in [0-100] before doing the fitting
#' @param verbose [boolean] Should warning messages about the data passed in be printed?
#' @return Plots to the active graphics device and returns and invisible NULL.
#' @export
#' @import RColorBrewer
#' @importFrom graphics plot rect
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom grDevices rgb
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @importFrom magicaxis magaxis


drugDoseResponseCurve <- 
function(drug, 
         cellline,
         pSets=list(),
         concentrations=list(),
         viabilities=list(), 
         conc_as_log = FALSE,
         viability_as_pct = TRUE,
         trunc=TRUE,
         legends.label = c("ic50_published", "gi50_published","auc_published","auc_recomputed","ic50_recomputed"),
         ylim=c(0,100), 
         xlim, mycol, 
         title,
         plot.type=c("Fitted","Actual", "Both"), 
         summarize.replicates=TRUE,
         lwd = 0.5,
         cex = 0.7,
         cex.main = 0.9, 
         legend.loc = "topright",
         verbose=TRUE) {
  if(!missing(pSets)){
    if (class(pSets) != "list") {
      if (class(pSets) == "PharmacoSet") {
        temp <- pSetName(pSets)
        pSets <- list(pSets)
        names(pSets) <- temp
      } else {
        stop("Type of pSets parameter should be either a pSet or a list of pSets.")
      }
    }
  }
  if(!missing(pSets) && (missing(drug) || missing(cellline))){
    stop("If you pass in a pSet then drug and cellline must be set") }
  # } else {
  #   if(missing(drug)){
  #   drug <- "Drug"}
  #   if(missing(cellline))
  #   cellline <- "Cell Line" 
  # }
  if(!missing(concentrations)){
    if(missing(viabilities)){

      stop("Please pass in the viabilities to Plot with the concentrations.")

    }
    if (class(concentrations) != "list") {
      if (mode(concentrations) == "numeric") {
        if(mode(viabilities)!="numeric"){
          stop("Passed in 1 vector of concentrations but the viabilities are not numeric!")
        }
        cleanData <- sanitizeInput(concentrations,
          viabilities,
          conc_as_log = conc_as_log,
          viability_as_pct = viability_as_pct,
          trunc = trunc,
          verbose = verbose)
        concentrations <- 10^cleanData[["log_conc"]]
        concentrations <- list(concentrations)
        viabilities <- 100*cleanData[["viability"]]
        viabilities <- list(viabilities)
        names(concentrations) <- "Exp1"
        names(viabilities) <- "Exp1"
      } else {
        stop("Mode of concentrations parameter should be either numeric or a list of numeric vectors")
      }
    } else{
      if(length(viabilities)!= length(concentrations)){
        stop("The number of concentration and viability vectors passed in differs")
      }
      if(is.null(names(concentrations))){
        names(concentrations) <- paste("Exp", seq_len(length(concentrations)))
      }
      for(i in seq_len(length(concentrations))){

        if (mode(concentrations[[i]]) == "numeric") {
          if(mode(viabilities[[i]])!="numeric"){
            stop(sprintf("concentrations[[%d]] are numeric but the viabilities[[%d]] are not numeric!",i,i))
          }
          cleanData <- sanitizeInput(concentrations[[i]],
            viabilities[[i]],
            conc_as_log = conc_as_log,
            viability_as_pct = viability_as_pct,
            trunc = trunc,
            verbose = verbose)
          concentrations[[i]] <- 10^cleanData[["log_conc"]]
          viabilities[[i]] <- 100*cleanData[["viability"]]
        } else {
          stop(sprintf("Mode of concentrations[[%d]] parameter should be numeric",i))
        }

      }

    }
  }

  common.range.star <- FALSE

  if (missing(plot.type)) {
    plot.type <- "Actual"
  }

  doses <- list(); responses <- list(); legend.values <- list(); j <- 0; pSetNames <- list()
  if(!missing(pSets)){
    for(i in seq_len(length(pSets))) {
      exp_i <- which(sensitivityInfo(pSets[[i]])[ ,"cellid"] == cellline & sensitivityInfo(pSets[[i]])[ ,"drugid"] == drug)
      if(length(exp_i) > 0) {
        if (summarize.replicates) {
          pSetNames[[i]] <- pSetName(pSets[[i]])
          if (length(exp_i) == 1) {
            drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp_i, , "Dose"])),
              "Viability"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp_i, , "Viability"])), stringsAsFactors=FALSE))
            drug.responses <- drug.responses[complete.cases(drug.responses), ]
          }else{
            drug.responses <- as.data.frame(cbind("Dose"=apply(pSets[[i]]@sensitivity$raw[exp_i, , "Dose"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}),
              "Viability"=apply(pSets[[i]]@sensitivity$raw[exp_i, , "Viability"], 2, function(x){median(as.numeric(x), na.rm=TRUE)}), stringsAsFactors=FALSE))
            drug.responses <- drug.responses[complete.cases(drug.responses), ]
          }
          doses[[i]] <- drug.responses$Dose
          responses[[i]] <- drug.responses$Viability
          names(doses[[i]]) <- names(responses[[i]]) <- seq_len(length(doses[[i]]))
          if (!missing(legends.label)) {
            if (length(legends.label) > 1) {
              legend.values[[i]] <- paste(unlist(lapply(legends.label, function(x){
                sprintf("%s = %s", x, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp_i,x]), digits=2))
              })), collapse = ", ")
            } else {
              legend.values[[i]] <- sprintf("%s = %s", legends.label, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp_i, legends.label]), digits=2))
            }
          } else {
            legend.values[[i]] <- ""
          }
        }else {
          for (exp in exp_i) {
            j <- j + 1
            pSetNames[[j]] <- pSetName(pSets[[i]])

            drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp, , "Dose"])),
              "Viability"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp, , "Viability"])), stringsAsFactors=FALSE))
            drug.responses <- drug.responses[complete.cases(drug.responses), ]
            doses[[j]] <- drug.responses$Dose
            responses[[j]] <- drug.responses$Viability
            names(doses[[j]]) <- names(responses[[j]]) <- seq_len(length(doses[[j]]))
            if (!missing(legends.label)) {
              if (length(legends.label) > 1) {
                legend.values[[j]] <- paste(unlist(lapply(legends.label, function(x){
                  sprintf("%s = %s", x, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp, x]), digits=2))
                })), collapse = ", ")
              } else {
                legend.values[[j]] <- sprintf(" Exp %s %s = %s", rownames(pSets[[i]]@sensitivity$info)[exp], legends.label, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp, legends.label]), digits=2))
              }
            } else {
              tt <- unlist(strsplit(rownames(pSets[[i]]@sensitivity$info)[exp], split="_"))
              if (tt[1] == "drugid") {
                legend.values[[j]] <- tt[2]
              }else{
                legend.values[[j]] <- rownames(pSets[[i]]@sensitivity$info)[exp]
              }
            }
          }
        }
      } else {
        warning("The cell line and drug combo were not tested together. Aborting function.")
        return()
      }
    }
  }

  if(!missing(concentrations)){
    doses2 <- list(); responses2 <- list(); legend.values2 <- list(); j <- 0; pSetNames2 <- list();
    for (i in seq_len(length(concentrations))){
      doses2[[i]] <- concentrations[[i]]
      responses2[[i]] <- viabilities[[i]]
      if(length(legends.label)>0){
        if(any(grepl("AUC", x=toupper(legends.label)))){
          legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "AUC", round(computeAUC(concentrations[[i]],viabilities[[i]], conc_as_log=FALSE, viability_as_pct=TRUE)/100, digits=2)), sep=", ")
        }
        if(any(grepl("IC50", x=toupper(legends.label)))){
          legend.values2[[i]] <- paste(legend.values2[i][[1]],sprintf("%s = %s", "IC50", round(computeIC50(concentrations[[i]],viabilities[[i]], conc_as_log=FALSE, viability_as_pct=TRUE), digits=2)), sep=", ")
        }

      } else{ legend.values2[[i]] <- ""}
      
      pSetNames2[[i]] <- names(concentrations)[[i]]
    }
    doses <- c(doses, doses2)
    responses <- c(responses, responses2)
    legend.values <- c(legend.values, legend.values2)
    pSetNames <- c(pSetNames, pSetNames2)
  }

  if (missing(mycol)) {
    # require(RColorBrewer) || stop("Library RColorBrewer is not available!")
    mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
  }

  dose.range <- c(10^100 , 0)
  viability.range <- c(0 , 10)
  for(i in seq_len(length(doses))) {
    dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
    viability.range <- c(0, max(viability.range[2], max(responses[[i]], na.rm=TRUE), na.rm=TRUE))
  }
  x1 <- 10 ^ 10; x2 <- 0

  if(length(doses) > 1) {
    common.ranges <- .getCommonConcentrationRange(doses)

    for(i in seq_len(length(doses))) {
      x1 <- min(x1, min(common.ranges[[i]]))
      x2 <- max(x2, max(common.ranges[[i]]))
    }
  }
  if (!missing(xlim)) {
    dose.range <- xlim
  }
  if (!missing(ylim)) {
    viability.range <- ylim
  }
  if(missing(title)){
    if(!missing(drug)&&!missing(cellline)){
      title <- sprintf("%s:%s", drug, cellline)
    } else {
      title <- "Drug Dose Response Curve"
    }
    
  }
  plot(NA, xlab="Concentration (uM)", ylab="% Viability", axes =FALSE, main=title, log="x", ylim=viability.range, xlim=dose.range, cex=cex, cex.main=cex.main)
  magicaxis::magaxis(side=seq_len(2), frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
  legends <- NULL
  legends.col <- NULL
  if (length(doses) > 1) {
    rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
  }

  for (i in seq_len(length(doses))) {
    points(doses[[i]],responses[[i]],pch=20,col = mycol[i], cex=cex)

    switch(plot.type , "Actual"={
      lines(doses[[i]], responses[[i]], lty=1, lwd=lwd, col=mycol[i])
    }, "Fitted"={ 
      log_logistic_params <- logLogisticRegression(conc=doses[[i]], viability=responses[[i]])
      log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
      lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=lwd, col=mycol[i])
    },"Both"={
      lines(doses[[i]],responses[[i]],lty=1,lwd=lwd,col = mycol[i])
      log_logistic_params <- logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
      log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
      lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=lwd, col=mycol[i])
    })
    legends<- c(legends, sprintf("%s%s", pSetNames[[i]], legend.values[[i]]))
    legends.col <-  c(legends.col, mycol[i])
  }
  if (common.range.star) {
    if (length(doses) > 1) {
      for (i in seq_len(length(doses))) {
        points(common.ranges[[i]], responses[[i]][names(common.ranges[[i]])], pch=8, col=mycol[i])
      }
    }
  }
  legend(legend.loc, legend=legends, col=legends.col, bty="n", cex=cex, pch=c(15,15))
  return(invisible(NULL))
}

