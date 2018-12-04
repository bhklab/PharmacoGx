
require(matrixStats)


merge_two_sides <- function(left, right, outx){
  
  left_vals <- left[[1]]
  left_discordant <- left[[2]]
  left_pairs <- left[[3]]
  
  right_vals <- right[[1]]
  right_discordant <- right[[2]]
  right_pairs <- right[[3]]
  
  RLR <- 0
  LLL <- length(left_vals)
  
  LR <- length(right_vals)
  
  out_vals <- numeric(LLL + LR)
  out_discordant <- numeric(length(out_vals))
  out_pairs <- numeric(length(out_vals))
  
  Li <- 1
  Ri <- 1
  i <- 1
  while(i <= length(out_vals)){
    
    if(LLL == 0){
      #Break out of loop if left list is empty
      out_vals[i] <- right_vals[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      Ri <- Ri + 1
      i <- i + 1
      next
    }
    if(RLR == LR){
      #Break out of loop if right list is empty
      out_vals[i] <- left_vals[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      Li <- Li + 1
      i <- i + 1
      next
    }
    
    if(left_vals[Li] < right_vals[Ri]) {
      out_vals[i] <- left_vals[Li]
      out_discordant[i] <- left_discordant[Li] + RLR
      out_pairs[i] <- left_pairs[Li]
      LLL <- LLL - 1
      Li <- Li + 1
      i <- i + 1
    } else if(left_vals[Li] > right_vals[Ri]) {
      out_vals[i] <- right_vals[Ri]
      out_discordant[i] <- right_discordant[Ri] + LLL
      out_pairs[i] <- right_pairs[Ri]
      RLR <- RLR + 1
      Ri <- Ri + 1
      i <- i + 1
    } else {
      # only case left is if the two values are equal.
      if(outx){
        out_vals[i] <- left_vals[Li]
        out_discordant[i] <- left_discordant[Li] + RLR 
        out_pairs[i] <- left_pairs[Li] - 1
        i <- i + 1
        out_vals[i] <- right_vals[Ri]
        out_discordant[i] <- right_discordant[Ri] + LLL - 1
        out_pairs[i] <- right_pairs[Ri] - 1
        LLL <- LLL - 1
        Li <- Li + 1
        RLR <- RLR + 1
        Ri <- Ri + 1
        i <- i + 1
      } else {
        # stop("Not implemented correctly?")
        out_vals[i] <- left_vals[Li]
        out_discordant[i] <- left_discordant[Li] + RLR + 0.5
        out_pairs[i] <- left_pairs[Li]
        i <- i + 1
        out_vals[i] <- right_vals[Ri]
        out_discordant[i] <- right_discordant[Ri] + LLL - 0.5
        out_pairs[i] <- right_pairs[Ri]
        LLL <- LLL - 1
        Li <- Li + 1
        RLR <- RLR + 1
        Ri <- Ri + 1
        i <- i + 1
      }
    }
  }
  
  return(list(out_vals, out_discordant, out_pairs))
   
}


merge_sort <- function(input, outx){
  if(length(input[[1]]) == 1){
    return(input)
  } else {
    input_vals <- input[[1]]
    input_discordant <- input[[2]]
    input_pairs <- input[[3]]
    split_idx <- floor(length(input_vals)/2)
    left <- list(input_vals[seq(1, split_idx)], 
                 input_discordant[seq(1, split_idx)],
                 input_pairs[seq(1, split_idx)])
    right <- list(input_vals[seq(split_idx+1, length(input_vals))], 
                  input_discordant[seq(split_idx+1, length(input_vals))],
                  input_pairs[seq(split_idx+1, length(input_vals))])
    left <- merge_sort(left, outx)
    right <- merge_sort(right, outx)
    output <- merge_two_sides(left, right, outx)
    return(output)
  }
}
## TODO: this code does not handle missing values well. 

fastCI <- function(observations, predictions, outx = TRUE, alpha = 0.05, alternative = c("two.sided", "greater", "less")){

  alternative = match.arg(alternative)
  if(!length(observations) == length(predictions)){
    stop("Size of vectors must be the same")
  }
  
  
  myCompleteCases <- complete.cases(observations, predictions)
  observations <- observations[myCompleteCases]
  predictions <- predictions[myCompleteCases]

  
  myorder <- order(observations)
  
  predictions <- predictions[myorder]

  input <- list(predictions, numeric(length(predictions)), rep(length(predictions)-1, length(predictions)))
  output <- merge_sort(input, outx)
  output_discordant <- output[[2]]
  output_pairs <- output[[3]]
  comppairs=10
  # N <- length(predictions)
  # D <- exp(logSumExp(log(output_discordant)))
  # C <- exp(logSumExp(log((N-1)-output_discordant)))
  # CC <- exp(logSumExp(log(C) + log(C-1)))
  # DD <- exp(logSumExp(log(D) + log(D-1)))
  # CD <- exp(logSumExp(log(C) + log(D)))

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
