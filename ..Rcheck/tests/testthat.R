Sys.unsetenv("R_TESTS")

library(testthat)
library(PharmacoGx)

test_check("PharmacoGx")
