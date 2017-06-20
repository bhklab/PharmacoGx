Sys.unsetenv("R_TESTS")

library(covr)
covr::codecov(quiet = FALSE)