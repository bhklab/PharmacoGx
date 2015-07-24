########################
## Zhaleh Safikhani
## All rights Reserved
## May 26, 2015
## Function to plot drug response curves for a list of pSets
########################


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
#' @param plot_type [character] Plot type which can be the actual one ("Actual") or 
#' the one fitted by logl logistic regression ("Fitted") or both of them ("Both").
#' If this parameter is missed by default actual curve is plotted. 
#' @import RColorBrewer
#' @import magicaxis
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom grDevices rgb
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @export


drugDoseResponseCurve <- 
  function(drug, cellline, pSets = list(), legends.label = c("ic50_published", "gi50_published","auc_published","auc_recomputed","ic50_recomputed"), ylim, xlim, mycol, plot_type=c("Fitted","Actual", "Both")) {
    if (typeof(pSets) != "list"){stop("type of pSets parameter should be list.")}
    #require(magicaxis) || stop("Library magicaxis is not available!")
    
    if(missing(plot_type)){
      plot_type <- "Actual"
    }
    
    doses <- list(); responses <- list(); legend.values <- list()
    for(i in 1:length(pSets))
    {
      exp_i <- which(pSets[[i]]@sensitivity$info[,"cellid"] == cellline & pSets[[i]]@sensitivity$info[,"drugid"] == drug)
      if(length(exp_i) > 0)
      {
        doses[[i]] <- as.numeric(pSets[[i]]@sensitivity$raw[exp_i, ,"Dose"])
        responses[[i]] <- as.numeric(pSets[[i]]@sensitivity$raw[exp_i, ,"Viability"])
        names(doses[[i]]) <- names(responses[[i]]) <- 1:length(doses[[i]])
        if(!missing(legends.label))
        {
          legend.values[[i]] <- paste(unlist(lapply(legends.label, function(x){
            sprintf("%s = %s", x, round(as.numeric(pSets[[i]]@sensitivity$phenotype[paste(drug, cellline, sep = "_"),x]), digits = 2))
          })), collapse = " , ")
        }
      }
    }
    
    if(missing(mycol))
    {
      #require(RColorBrewer) || stop("Library RColorBrewer is not available!")
      mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
    }
    
    dose.range <- c(10^100 , 0)
    viability.range <- c(0 , 10)
    for(i in 1:length(doses))
    {
      dose.range <- c(min(dose.range[1],min(doses[[i]], na.rm = T), na.rm = T), max(dose.range[2],max(doses[[i]], na.rm = T), na.rm = T))
      viability.range <- c(0, max(viability.range[2],max(responses[[i]], na.rm = T), na.rm = T))
    }
    x1 <- 10^10; x2 <- 0
    
    if(length(doses) > 1)
    {
      common.ranges <-  .getCommonConcentrationRange(doses)
      
      for(i in 1:length(doses))
      {
        x1 <- min(x1, min(common.ranges[[i]]))
        x2 <- max(x2, max(common.ranges[[i]]))
      }
    }
    if(!missing(xlim)){dose.range <- xlim}
    if(!missing(ylim)){viability.range <- ylim}
    
    plot(NA, xlab = "Concentration (uM)",ylab="% Viability",axes =FALSE, main=sprintf("%s:%s",drug,cellline),log="x",ylim = viability.range, xlim = dose.range, pch=8,cex=.7, cex.main = .9)
    magicaxis::magaxis(side=1:2,frame.plot = TRUE,tcl=-.3,majorn=c(5,3),minorn=c(5,2))
    legends <- NULL
    legends.col <- NULL
    if(length(doses) > 1)
    {
      rect(xleft = x1  , xright = x2 , ybottom = viability.range[1] , ytop = viability.range[2] , col = rgb(240,240,240, maxColorValue = 255), border = F)
    }
    
    for(i in 1:length(doses))
    {
      points(doses[[i]],responses[[i]],pch=20,col = mycol[i])
      
      switch(plot_type , "Actual" = {
        lines(doses[[i]],responses[[i]],lty=1,lwd=.5,col = mycol[i])
        }, "Fitted" = { 
          log_logistic_params <- logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals,.Hill(log10_x_vals, pars = c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1,lwd=.5,col = mycol[i])
        },"Both"={
          lines(doses[[i]],responses[[i]],lty=1,lwd=.5,col = mycol[i])
          log_logistic_params <- logLogisticRegression(conc = doses[[i]], viability = responses[[i]])
          log10_x_vals <- .GetSupportVec(log10(doses[[i]]))
          lines(10 ^ log10_x_vals,.Hill(log10_x_vals, pars = c(log_logistic_params$HS, log_logistic_params$E_inf/100, log10(log_logistic_params$EC50))) * 100 ,lty=1,lwd=.5,col = mycol[i])
        })
      if(!missing(legends.label))
      {
        legends<- c(legends,sprintf("%s %s", pSetName(pSets[[i]]), legend.values[[i]]))
        legends.col <-  c(legends.col, mycol[i])
      }
    }
    if(length(doses) > 1)
    {
      for(i in 1:length(doses))
      {
        points(common.ranges[[i]],responses[[i]][names(common.ranges[[i]])],pch=8,col = mycol[i])
      }
    }
    
    if(!missing(legends.label))
    {
      legend("topright",legend=legends, col=legends.col, bty="n", cex = .7, pch = c(15,15))
    }
    
  }

