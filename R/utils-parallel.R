# Some restricted runtimes cannot detect a CPU count. Keep signature
# evaluation usable and avoid starting workers without a verified limit.
.capSignatureThreads <- function(nthread, available = parallel::detectCores()) {
  if (is.na(available)) {
    available <- 1L
  }
  if (nthread > available) available else nthread
}
