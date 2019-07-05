context("Homozygous deletion and duplication")

test_that("summarize CNP_029", {
  library(panc.data)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  if(!dir.exists(path)) dir.create(path)
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[29, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
})

test_that("Homozygous deletion / duplication pipeline", {
  library(SummarizedExperiment)
  set.seed(5)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_029"]
  snpdat <- snp_se[overlapsAny(snp_se, g), ]e
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  THR <- summaries(mb.subsamp)$deletion_cutoff
  simdat <- augment_homozygous(mb.subsamp)
  sb <- warmup(assays(mb.subsamp),
               "SBP4",
               "SB4")
  mp <- McmcParams(iter=400, burnin=500)
  mcmcParams(sb) <- mp
  sb <- posteriorSimulation(sb)
  finished <- stop_early(sb, 0.99, 0.99)
  expect_false(finished)
  THR <- summaries(mb.subsamp)$deletion_cutoff
  ##trace(explore_multibatch, browser)
  mb <- revertToMultiBatch(sb)
  restricted <- fit_restricted(mb, sb, THR, model="MBP3")
  expect_true(any(round(p(restricted), 2)[, 3] < 0.05))
  ## duplications are rare and will cause problems
  dat <- filter(assays(sb), oned > 0.1)
  cutoff <- duplication_cutoff(assays(sb))
  ##
  ## impute rare duplications
  ##
  ## fit model
  ##
  ## Need unit tests on augmentation
  ## rare homozygous
  ## rare hemizygous
  ## rare duplications
  ## rare hemizygous and rare duplications

  mb <- explore_multibatch(sb, simdat, THR,
                           model="MBP3")
  choose_mb <- log_lik(mb) > log_lik(sb) + 25
  if(choose_mb) full <- mb else  full <- sb
  gmodel <- genotype_model(full, snpdat)

})
