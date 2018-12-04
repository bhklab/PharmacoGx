source("~/Code/fastCI.R")

fastCIByGroup <- function(observations, predictions, groups, outx = TRUE, alpha = 0.05, alternative = c("two.sided", "greater", "less")){
  
  
  alternative = match.arg(alternative)
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  if(!unique(c(length(observations), length(predictions), length(groups))) ==  length(observations)){
    stop("Please pass in vectors of all the same length")
  }
  
  myCompleteCases <- complete.cases(observations, predictions, groups)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]
  groups <- groups[myCompleteCases]
  
  myobs <- split(observations, as.factor(groups))
  mypreds <- split(predictions, as.factor(groups))
  
  
  res <- mapply(function(observations, predictions, outx, alternative){
    # browser()
    myorder <- order(observations)
    predictions <- predictions[myorder]
    input <- list(predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
    output <- merge_sort(input, outx)
    return(output)
  },  myobs, mypreds, outx)

  output_discordant <- unlist(res[2,])
  output_pairs <- unlist(res[3,])
  comppairs=10
  
  N <- length(predictions)
  D <- sum(output_discordant)
  Cvec <- output_pairs-output_discordant
  C <-  sum(Cvec)
  CC <- sum(Cvec*(Cvec-1))
  DD <- sum(output_discordant*(output_discordant-1))
  CD <- sum(Cvec*output_discordant)
  # browser()
  if (N < 3 || (C == 0 && D == 0)) {
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=0))
  }
  if(C==0 || D==0 || C * (C - 1)==0 || D * (D - 1)==0 || C * D==0 || (C + D) < comppairs){
    return(list("cindex"=NA, "p.value"=NA, "sterr"=NA, "lower"=NA, "upper"=NA, "relevant.pairs.no"=(C + D) / 2))
  }
  # cindex <- exp(C) / exp(logSumExp(c(C, D)))
  cindex <- C/(C+D)
  varp <- 4 * ((D ^ 2 * CC - 2 * C * D * CD + C ^ 2 * DD) / (C + D) ^ 4) * N * (N - 1) / (N - 2) 
  
  # varp <- 4 * ((exp(logSumExp(c(2*D + CC, 2*C + DD))) - 2 *exp(C + D + CD)) / exp(logSumExp(c(C, D)))^4) * N * (N - 1) / (N - 2)
  
  if (varp >= 0) {
    sterr <- sqrt(varp / N)
    ci <- qnorm(p = alpha / 2, lower.tail = FALSE) * sterr
    p <- pnorm((cindex - 0.5) / sterr)
  } else {
    return(list("cindex"=cindex,
                "p.value"=1,
                "sterr"=NA,
                "lower"=0,
                "upper"=0,
                "relevant.pairs.no"=(C + D) / 2))
  }
  return(list("cindex"=cindex,
              "p.value"=switch(alternative, less=p, greater=1 - p, two.sided=2 * min(p, 1 - p)),
              "sterr"=sterr,
              "lower"=max(cindex - ci, 0),
              "upper"=min(cindex + ci, 1),
              "relevant.pairs.no"=(C + D) / 2))
}