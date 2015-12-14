#' Plot drug response curve of a given drug and a given cell for a list of pSets (objects of the PharmacoSet class).
#' 
#' Given a list of PharmacoSets, the function will plot the drug_response curve,
#' for a given drug/cell pair. The y axis of the plot is the viability percentage
#' and x axis is the log transformed concentrations. If more than one pSet is 
#' provided, a light gray area would show the common concentration range between pSets. 
#' User can ask for type of sensitivity measurment to be shown in the plot legend.
#' PharmacoSet object's internal curation slot, or according to matching tables 
#' provided by the user
#' 
#' @param drug [string] A drug name for which the drug response curve should be 
#' plotted. If the plot is desirable for more than one pharmaco set, A unique drug id
#' should be provided.
#' @param cellline [string] A cell line name for which the drug response curve should be 
#' plotted. If the plot is desirable for more than one pharmaco set, A unique cell id
#' should be provided.
#' @param pSets [list] a list of PharmacoSet objects, for which the function
#' should plot the curves.
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
  function(drug, cellline, pSets=list(), legends.label = c("ic50_published", "gi50_published","auc_published","auc_recomputed","ic50_recomputed"), ylim, xlim, mycol, plot.type=c("Fitted","Actual", "Both"), summarize.replicates=TRUE) {
    
    if (class(pSets) != "list") {
      if (class(pSets) == "PharmacoSet") {
        temp <- pSetName(pSets)
        pSets <- list(pSets)
        names(pSets) <- temp
      } else {
        stop("type of pSets parameter should be either a pSet or a list of pSets.")
      }
    }
    common.range.star <- FALSE
    
    if (missing(plot.type)) {
      plot.type <- "Actual"
    }
    
    doses <- list(); responses <- list(); legend.values <- list(); j <- 0; pSetIndex <- list()
    for(i in 1:length(pSets)) {
      exp_i <- which(sensitivityInfo(pSets[[i]])[ ,"cellid"] == cellline & sensitivityInfo(pSets[[i]])[ ,"drugid"] == drug)
      if(length(exp_i) > 0) {
        if (summarize.replicates) {
          pSetIndex[[i]] <- i
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
          names(doses[[i]]) <- names(responses[[i]]) <- 1:length(doses[[i]])
          if (!missing(legends.label)) {
            if (length(legends.label) > 1) {
              legend.values[[i]] <- paste(unlist(lapply(legends.label, function(x){
                sprintf("%s = %s", x, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp_i,x]), digits=2))
                })), collapse = " , ")
            } else {
              legend.values[[i]] <- sprintf("%s = %s", legends.label, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp_i, legends.label]), digits=2))
            }
          } else {
            legend.values[[i]] <- ""
          }
        }else {
          for (exp in exp_i) {
            j <- j + 1
            pSetIndex[[j]] <- i
            
            drug.responses <- as.data.frame(cbind("Dose"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp, , "Dose"])),
                                                  "Viability"=as.numeric(as.vector(pSets[[i]]@sensitivity$raw[exp, , "Viability"])), stringsAsFactors=FALSE))
            drug.responses <- drug.responses[complete.cases(drug.responses), ]
            doses[[j]] <- drug.responses$Dose
            responses[[j]] <- drug.responses$Viability
            names(doses[[j]]) <- names(responses[[j]]) <- 1:length(doses[[j]])
            if (!missing(legends.label)) {
              if (length(legends.label) > 1) {
                legend.values[[j]] <- paste(unlist(lapply(legends.label, function(x){
                  sprintf("%s = %s", x, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp, x]), digits=2))
                })), collapse = " , ")
              } else {
                legend.values[[j]] <- sprintf("Exp %s %s = %s", rownames(pSets[[i]]@sensitivity$info)[exp], legends.label, round(as.numeric(pSets[[i]]@sensitivity$profiles[exp, legends.label]), digits=2))
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
    }
  }
    
    if (missing(mycol)) {
#       require(RColorBrewer) || stop("Library RColorBrewer is not available!")
      mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
    }
    
    dose.range <- c(10^100 , 0)
    viability.range <- c(0 , 10)
    for(i in 1:length(doses)) {
      dose.range <- c(min(dose.range[1], min(doses[[i]], na.rm=TRUE), na.rm=TRUE), max(dose.range[2], max(doses[[i]], na.rm=TRUE), na.rm=TRUE))
      viability.range <- c(0, max(viability.range[2], max(responses[[i]], na.rm=TRUE), na.rm=TRUE))
    }
    x1 <- 10 ^ 10; x2 <- 0
    
    if(length(doses) > 1) {
      common.ranges <- .getCommonConcentrationRange(doses)
      
      for(i in 1:length(doses)) {
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
    
    plot(NA, xlab="Concentration (uM)", ylab="% Viability", axes =FALSE, main=sprintf("%s:%s", drug, cellline), log="x", ylim=viability.range, xlim=dose.range, cex=.7, cex.main=.9)
    magicaxis::magaxis(side=1:2, frame.plot=TRUE, tcl=-.3, majorn=c(5,3), minorn=c(5,2))
    legends <- NULL
    legends.col <- NULL
    if (length(doses) > 1) {
      rect(xleft=x1, xright=x2, ybottom=viability.range[1] , ytop=viability.range[2] , col=rgb(240, 240, 240, maxColorValue = 255), border=FALSE)
    }
    
    for (i in 1:length(doses)) {
      points(doses[[i]],responses[[i]],pch=20,col = mycol[i])
      
      switch(plot.type , "Actual"={
        lines(doses[[i]], responses[[i]], lty=1, lwd=.5, col=mycol[i])
        }, "Fitted"={ 
          log_logistic_params <- logLogisticRegression(conc=doses[[i]], viability=responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=.5, col=mycol[i])
        },"Both"={
          lines(doses[[i]],responses[[i]],lty=1,lwd=.5,col = mycol[i])
          log_logistic_params <- logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals, .Hill(log10_x_vals, pars=c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1, lwd=.5, col=mycol[i])
        })
        legends<- c(legends, sprintf("%s %s", pSetName(pSets[[pSetIndex[[i]]]]), legend.values[[i]]))
        legends.col <-  c(legends.col, mycol[i])
    }
    if (common.range.star) {
      if (length(doses) > 1) {
        for (i in 1:length(doses)) {
          points(common.ranges[[i]], responses[[i]][names(common.ranges[[i]])], pch=8, col=mycol[i])
        }
      }
    }
      legend("topright", legend=legends, col=legends.col, bty="n", cex=.7, pch=c(15,15))
  }

