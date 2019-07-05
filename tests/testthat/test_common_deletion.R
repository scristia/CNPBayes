context("Common deletions")

test_that("summarize_cnp22", {
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_022")
  set.seed(555)
  library(panc.data)
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[22, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
})

test_that("common deletion", {
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_022")
  CNP <- "CNP_022"
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  sb3 <- warmup(assays(mb.subsamp),
               "SBP3",
               "SB3")
  mp <- McmcParams(iter=400, burnin=500)
  mcmcParams(sb3) <- mp
  sb3 <- posteriorSimulation(sb3)
  finished <- stop_early(sb3)
  expect_false(finished)
  mb <- warmup(assays(mb.subsamp),
               "MBP3", "MB3")
  mcmcParams(mb) <- mp
  mb <- posteriorSimulation(mb)
  finished <- stop_early(mb)
  expect_false(finished)
  mod_1.3 <- explore_restricted(mb, sb3, -1)
  if(FALSE){
    saveRDS(mod_1.3, file=file.path(path, paste0(CNP, ".rds")))
  }
  expected <- readRDS(file.path(path, paste0(CNP, ".rds")))
  expect_equivalent(mod_1.3, expected)
})
