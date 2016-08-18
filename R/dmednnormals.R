
.multinom<-function(x, y) {
  coeff <- 1
  for (i in 1:length(y)) {
    coeff <- coeff * choose(x, y[i])
    x <- x - y[i]
  }
  return(coeff)
}

.rmednnormals = function(N, n, scale) {
  x <- matrix(NA, nrow = 1, ncol = N)
  for (i in 1:N) {
    x[i] <- median(rnorm(n, sd = scale))
  }
  return(x)
}

.dmednnormals = function(x, n, scale, divisions = 100) {
  n <- rep(n, times = length(x) / length(n))
  scale <- rep(scale, times = length(x) / length(scale))
  y <- matrix(NA, nrow = 1, ncol = length(x))
  for (g in seq_along(x)) {
    if (n[g] %% 2 == 0) {
      y[g] <- 2 * .multinom(n[g], c(n[g] / 2 - 1, n[g] / 2 - 1)) *
        tryCatch(integrate(f = function(j) {
          (pnorm(x[g] - j / 2, sd = scale[g])) ^ (n[g] / 2 - 1) *
            (1 - pnorm(x[g] + j / 2, sd = scale[g])) ^ (n[g] / 2 - 1) *
            dnorm(x[g] - j / 2, sd = scale[g]) *
            dnorm(x[g] + j / 2, sd = scale[g])},
          lower = 0,
          upper = Inf,
          subdivisions = divisions)[[1]],
          error = function(e) {
            if (divisions == 1) {
              wseq <- c(1, 4, 1)
            } else {
              wseq <- c(1, 4, rep(c(2, 4), times = divisions - 1), 1)
            }
            aseq <- seq(from = 0, to = pi / 2, length.out = 2 * divisions + 1)
            tseq <- tan(aseq) / 2
            return(sum((pnorm(x[g] + tseq, sd = scale[g])) ^ (n[g] / 2 - 1) *
                         (pnorm(x[g] - tseq, sd = scale[g])) ^ (n[g] / 2 - 1) *
                         dnorm(x[g] + tseq, sd = scale[g]) * 
                         dnorm(x[g] - tseq, sd = scale[g]) /
                         (cos(aseq)) ^ 2 * wseq) *
                     (aseq[2] - aseq[1]) / 6)
          })
    } else {
      y[g] <- .multinom(n[g], c((n[g] - 1) / 2, (n[g] - 1) / 2)) *
        (pnorm(x[g], sd = scale[g])) ^ ((n[g] - 1) / 2) *
        (1 - pnorm(x[g], sd = scale[g])) ^ ((n[g] - 1) / 2) *
        dnorm(x[g], sd = scale[g])
    }
  }
  return(y)
}

.pmednnormals = function(x, n, scale, divisions = 100) {
  n <- rep(n, times = length(x) / length(n))
  scale <- rep(scale, times = length(x) / length(scale))
  y <- numeric(length(x))
  for (g in seq_along(x)) {
    if (n[g] %% 2 == 0) {
      y[g] <- tryCatch(integrate(f = function(k) {.dmednnormals(k, n[g], scale[g])},
                                 lower = -Inf,
                                 upper = x[g],
                                 subdivisions = divisions)[[1]],
                       error = function(e) {
                         wseq <- c(1, 4, rep(c(2, 4), times = divisions - 1), 1)
                         aseq <- seq(from = -pi / 2, to = atan(x[g]), length.out = 2 * divisions + 1)
                         return(sum(.dmednnormals(tan(aseq), n[g], scale[g]) * wseq / (cos(aseq)) ^ 2) * (aseq[3] - aseq[1]) / 6)})
    } else {
      y[g] <- 0
      Fx <- pnorm(x[g], sd = scale[g])
      for (i in 0:((n[g] - 1) / 2)) {
        y[g] <- y[g] + choose((n[g] - 1) / 2, i) * (-1) ^ i * Fx ^ ((n[g] + 1) / 2 + i) / ((n[g] + 1) / 2 + i)
      }
      y[g] <- y[g] * .multinom(n[g], c((n[g] - 1) / 2, (n[g] - 1) / 2))
    }
  }
  return(y)
}

.edmednnormals = function(x, n, scale, divisions = 100) {
  n <- rep(n, times = length(x) / length(n))
  scale <- rep(scale, times = length(x) / length(scale))
  y <- numeric(length(x))
  for (g in seq_along(y)) {
    if (x[g] > 0) {
      upper <- Inf
      lower <- x[g]
    } else {
      upper <- x[g]
      lower <- -Inf
    }
    y[g] <- tryCatch(integrate(f = function(k) {(.dmednnormals(k, n[g], scale[g])) ^ 2},
                               lower = lower,
                               upper = upper,
                               subdivisions = divisions)[[1]],
                     error = function(e) {
                       wseq <- c(1, 4, rep(c(2, 4), times = divisions - 1), 1)
                       aseq <- seq(from = atan(lower), to = atan(upper), length.out = 2 * divisions + 1)
                       return(sum((.dmednnormals(tan(aseq), n[g], scale[g])) ^ 2 * wseq / (cos(aseq)) ^ 2) * (aseq[3] - aseq[1]) / 6)})
  }
  return(y)
}