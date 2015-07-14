#HIDDEN FUNCTIONS

#PREDICT VIABILITY FROM CONCENTRATION DATA AND CURVE PARAMETERS
.Hill<-function(x, pars) {
  return(pars[2] + (1 - pars[2]) / (1 + (10 ^ x / 10 ^ pars[3]) ^ pars[1]))
}

#CALCULATE RESIDUAL OF FIT
.Residual<-function(x, y, pars, scale = 0.07, Cauchy_flag = TRUE, trunc = FALSE) {
  if (Cauchy_flag == FALSE) {
    return(sum((.Hill(x, pars) - y) ^ 2))
  } else {
    diffs <- .Hill(x, pars)-y
    if (trunc == FALSE) {
      return(sum(-log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs / scale)) * (1 / 2 - 1 / pi * atan(diffs / scale)))))
    } else {
      down_truncated <- which(abs(y) > 1)
      up_truncated <- which(abs(y) < 0)
      return(sum(log(6 * scale / (pi * (scale ^ 2 + diffs ^ 2)) * (1 / 2 + 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))
                     * (1 / 2 - 1 / pi * atan(diffs[setdiff(1:length(y), union(down_truncated, up_truncated))] / scale))),
                 -log(1 / 2 - 3 / (2 * pi) * atan((1 - diffs[down_truncated] - y[down_truncated]) / scale) + 2 / pi ^ 3 * (atan((1 - diffs[down_truncated] - y[down_truncated]) / scale)) ^ 3),
                 -log(-1 / 2 + 3 / (2 * pi) * atan((-diffs[up_truncated] - y[up_truncated]) / scale) - 2 / pi ^ 3 * (atan((- diffs[up_truncated] - y[up_truncated]) / scale)) ^ 3)))
    }
  }
}

#GENERATE AN INITIAL GUESS FOR DOSE-RESPONSE CURVE PARAMETERS BY EVALUATING THE RESIDUALS AT DIFFERENT LATTICE POINTS OF THE SEARCH SPACE
.MeshEval<-function(log_conc,
                    viability,
                    lower_bounds = c(0, 0, -6),
                    upper_bounds = c(4, 1, 6),
                    density = c(2, 10, 2),
                    scale = 0.07,
                    Cauchy_flag = TRUE,
                    trunc = FALSE) {
  guess_residual<-Inf
  for (i in seq(from = lower_bounds[1], to = upper_bounds[1], by = 1 / density[1])) {
    for (j in seq(from = lower_bounds[2], to = upper_bounds[2], by = 1 / density[2])) {
      for (k in seq(from = lower_bounds[3], to = upper_bounds[3], by = 1 / density[3])) {
        test_guess_residual <- .Residual(log_conc,
                                         viability,
                                         pars = c(i, j, k),
                                         scale = scale,
                                         Cauchy_flag = Cauchy_flag,
                                         trunc = trunc)
        if (test_guess_residual < guess_residual) {
          guess <- c(i, j, k)
          guess_residual <- test_guess_residual
        }
      }
    }
  }
  return(guess)
}