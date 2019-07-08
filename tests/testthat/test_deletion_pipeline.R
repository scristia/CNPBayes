context("analysis of rare deletions")

test_that("rare_deletion_pipeline", {
  path <- system.file("extdata", package="CNPBayes")
  set.seed(2463)
  snp_se <- readRDS(file.path(path, "snp_se.rds"))
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  simdat <- augment_homozygous(mb.subsamp)
  THR <- summaries(mb.subsamp)$deletion_cutoff
  sb3 <- warmup(simdat, "SBP3", "SB3")
  mcmcParams(sb3) <- McmcParams(iter=400, burnin=500)
  sb3 <- posteriorSimulation(sb3)
  expect_equal(as.numeric(theta(sb3)),
               c(-3.76, -0.22, 0.002),
               tolerance=0.01)
  finished <- stop_early(sb3)
  expect_false(finished)
  ## Since not finished, keep going
  ##trace(explore_multibatch, browser)
  mod_1.3 <- explore_multibatch(sb3, simdat, THR)
  gmodel <- genotype_model(mod_1.3, snp_se)
  if(FALSE){
    saveRDS(sb3, file=file.path(path, "sb3_1.rds"))
  }
  expected_sb3 <- readRDS(file.path(path,
                                    "sb3_1.rds"))
  expected_mod1.3 <- readRDS(file.path(path, "mod_1.3.rds"))
  expected_gt <- readRDS(file.path(path, "CNP_001.rds"))
  expect_equivalent(sb3, expected_sb3)
  expect_equivalent(mod_1.3, expected_mod1.3)
  expect_equivalent(gmodel, expected_gt)
})

## for first step, see test_summarized_region.R
.test_that <- function(nm, expr) NULL

test_that("augment data", {
  set.seed(2463)
  ##mp <- McmcParams(iter=500, burnin=400)
  path <- system.file("extdata", package="CNPBayes")
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  simdat <- augment_homozygous(mb.subsamp)
  if(FALSE){
    saveRDS(simdat, file=file.path(path, "simdat.rds"))
  }
  expected <- readRDS(file.path(path, "simdat.rds"))
  expected$homozygousdel_mean <- simdat$homozygousdel_mean
  expect_equivalent(simdat, expected)
})

test_that("Augment hemizygous", {
  path <- system.file("extdata", package="CNPBayes")
  ##
  ## Reminder that we only do this section if augment_hemizygous
  ## is TRUE
  ##
  mod_2.3 <- readRDS(file.path(path, "mod_2.3.rds"))
  sb3 <- readRDS(file.path(path, "sb3_1.rds"))
  mb <- augment_rarecomponent(mod_2.3, sb3)
  ##
  message("Run additional MCMC simulations on the augmented data")
  ##
  mod_2.3 <- posteriorSimulation(mb)
  if(FALSE){
    saveRDS(mod_2.3, file=file.path(path, "mod_2.32.rds"))
  }
  expected <- readRDS(file.path(path, "mod_2.32.rds"))
  expect_equivalent(mod_2.3, expected)
  expect_equivalent(assays(mod_2.3), assays(expected))
})

test_that("Fit multi-batch model with all components", {
  set.seed(2463)
  path <- system.file("extdata", package="CNPBayes")
  simdat <- readRDS(file.path(path, "simdat_2.rds"))
  mod_2.3 <- readRDS(file.path(path, "mod_2.3.rds"))
  sb3 <- readRDS(file.path(path, "sb3_1.rds"))
  mb <- revertToMultiBatch(sb3)
  fdat <- filter(assays(mb), oned > -1)  %>%
     mutate(is_simulated=FALSE)
  mb <- warmup(fdat, "MBP2")
  mcmcParams(mb) <- McmcParams(iter=500, burnin=50)
  mod_1.3 <- mcmcWithHomDel(mb, sb=sb3,
                            restricted_model=mod_2.3,
                            THR=-1)
  ##
  ##
  ##
  ok <- ok_model(mod_1.3, restricted_model=mod_2.3)
  expect_true(ok)
  if(FALSE){
    saveRDS(mod_1.3, file=file.path(path, "mod_1.3.rds"))
  }
  expected <- readRDS(file.path(path, "mod_1.3.rds"))
  expect_equivalent(mod_1.3, expected)
})

test_that("only update homozygous component", {
  path <- system.file("extdata", package="CNPBayes")
  set.seed(2463)
  simdat <- readRDS(file.path(path, "simdat_1.rds"))
  mod_1.3 <- readRDS(file.path(path, "mod_1.3.rds"))
  mod_2.3 <- readRDS(file.path(path, "mod_2.3.rds"))
  sb3 <- readRDS(file.path(path, "sb3_1.rds"))
  if(FALSE){
    saveRDS(mod_1.3, file=file.path(path, "mod_1.3_2.rds"))
  }
  model <- gsub("P", "", modelName(mod_1.3))
  mod_1.3 <- mcmcHomDelOnly(simdat, mod_2.3, sb3, model)
  expected <- readRDS(file.path(path, "mod_1.3_2.rds"))
  expect_equivalent(mod_1.3, expected)
})

test_that("genotype", {
  path <- system.file("extdata", package="CNPBayes")
  ##
  ## Genotype using high confidence subjects
  ##
  path <- system.file("extdata", package="CNPBayes")
  snp_se <- readRDS(file.path(path, "snp_se.rds"))
  mod_1.3 <- readRDS(file.path(path, "mod_1.3.rds"))
  mod_1.3 <- genotype_model(mod_1.3, snp_se)
  ##
  ##
  ##
  message("Write mapping to model file")
  expect_identical(mapping(mod_1.3), c("0", "2", "2"))
  if(FALSE){
    saveRDS(mod_1.3, file=file.path(path, "CNP_001.rds"))
  }
  expected <- readRDS(file.path(path, "CNP_001.rds"))
  expect_equivalent(mod_1.3, expected)
  th.expected <- theta(expected)
  expect_equal(theta(mod_1.3), th.expected, tolerance=0.05)
})

## Plotting code
.test_that("plots", {
  path <- system.file("extdata", package="CNPBayes")
  mod_1.3 <- readRDS(file.path(path, "CNP_001.rds"))
  snpdat <- readRDS(file.path(path, "snp_se.rds"))
  bafdat2 <- join_baf_oned(mod_1.3, snpdat)
  figs <- list_mixture_plots(mod_1.3, bafdat2)
  pdf(tempfile(), width=14, height=8)
  mixture_layout(figs, augmented=TRUE)
  dev.off()

  pdf(tempfile(), width=14, height=8)
  mixture_layout(figs, augmented=FALSE)
  dev.off()
})
