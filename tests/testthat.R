Sys.setenv("R_TESTS" = "")

library(testthat)
library(bayesboot)

test_check("bayesboot")
