#' Function computing connectivity scores between two signatures
#'
#' @inherit CoreGx::connectivityScore
#' @inheritParams CoreGx::connectivityScore
#'
#' @export
connectivityScore <-
  function(
    x,
    y,
    method = c("gsea", "fgsea", "gwc"),
    nperm = 1e4,
    nthread = 1,
    gwc.method = c("spearman", "pearson"),
    ...
  ) {
    method <- match.arg(method)
    gwc.method <- match.arg(gwc.method)
    CoreGx::connectivityScore(x, y, method, nperm, nthread, gwc.method, ...)
  }
