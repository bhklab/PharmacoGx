Sys.unsetenv("R_TESTS")

library(covr)
options(covr.fix_parallel_mcexit=TRUE)
covr::codecov(quiet = FALSE)