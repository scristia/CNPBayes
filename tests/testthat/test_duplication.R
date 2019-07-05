context("Duplications")

test_that("summarize_cnp244", {
  library(panc.data)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_244")
  if(!dir.exists(path)) dir.create(path)
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[244, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
})

test_that("duplication", {
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_244")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_244"]
  snpdat <- snp_se[overlapsAny(snp_se, g), ]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  sb <- warmup(assays(mb.subsamp),
               "SBP2",
               "SB2")
  mp <- McmcParams(iter=1000, burnin=500)
  mcmcParams(sb) <- mp
  sb <- posteriorSimulation(sb)
  finished <- stop_early(sb, 0.98, 0.98)
  expect_false(finished)
  ##
  ## Try MultiBatch
  ##
  mb <- warmup(assays(mb.subsamp), "MBP2")
  mcmcParams(mb) <- mp
  mb <- posteriorSimulation(mb)
  choose_mb <- log_lik(mb) > log_lik(sb) + 25
  if(choose_mb) full <- mb else  full <- sb
  gmodel <- genotype_model(full, snpdat)
  expect_identical(mapping(gmodel), c("2", "3"))
  if(FALSE){
    bafdat <- join_baf_oned(full, snpdat)
    figs <- list_mixture_plots(full, bafdat)
    mixture_layout(figs, augmented=TRUE)
  }
})
