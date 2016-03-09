context("Bioconductor Check")

test_that("BiocCheck", {
    skip("Not running BiocCheck right now")
    # run BiocCheck on package
    library(BiocCheck)
    curr <- getwd()
    setwd("../..")
    check <- BiocCheck(".")

    # pull out the errors/warnings/messages into separate variables
    requirements <- check$requirements
    recommendations <- check$recommendations
    considerations <- check$considerations

    # warnings that are known and ignored
    ignore <- c("Remove vignette sources from inst/doc; they belong in vignettes/.",
                "y of x.y.z version should be even in release",
                "Register native routines! see http://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines")
    recommendations <- recommendations[!recommendations %in% ignore]

    # check that there are no errors
    expect_true(length(requirements) == 0)
    expect_true(length(recommendations) == 0)
    setwd(curr)
})
