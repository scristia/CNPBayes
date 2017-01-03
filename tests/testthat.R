Sys.setenv("R_TESTS"="")
library(testthat)
library(CNPBayes)
test_check("CNPBayes")
