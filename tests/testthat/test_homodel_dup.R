context("Homozygous deletion and duplication")

.test_that <- function(nm, expr) NULL

test_that("summarize CNP_029", {
  set.seed(5)
  library(panc.data)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  if(!dir.exists(path)) dir.create(path)
  data(cnp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[29, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
  expect_equivalent(assays(mb.subsamp), assays(expected))
})

test_that("Homozygous deletion / duplication pipeline", {
  library(SummarizedExperiment)
  set.seed(5)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_029"]
  snpdat <- snp_se[overlapsAny(snp_se, g), ]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  THR <- summaries(mb.subsamp)$deletion_cutoff
  simdat <- augment_homozygous(mb.subsamp)
  sb <- warmup(assays(mb.subsamp),
               "SBP4",
               "SB4")
  mp <- McmcParams(iter=400, burnin=500)
  mcmcParams(sb) <- mp
  sb <- posteriorSimulation(sb)
  ## at least 99% of samples above 0.95 probability
  finished <- stop_early(sb, 0.95, 0.99)
  expect_false(finished)
  expect_true(equivalent_variance(sb))
  THR <- summaries(mb.subsamp)$deletion_cutoff
  ##
  ## even though we could have stopped, try to improve by multibatch
  ##
  fdat <- filter(assays(mb.subsamp), oned > THR)
  mb <- warmup(fdat, "MBP3")
  mcmcParams(mb) <- mp
  set.seed(123)
  mod_2.4 <- restricted_homhemdup(mb, mod_2.4, mb.subsamp, mp)
  if(FALSE){
    ##fig1=ggMixture(tmp)
    ##fig2=ggMixture(mod_2.4)
    grid.arrange(fig1, fig2)
    saveRDS(mod_2.4, file=file.path(path, "mod_2.4.rds"))
  }
  expected <- readRDS(file.path(path, "mod_2.4.rds"))
  expect_equivalent(assays(mod_2.4), assays(expected))
  expect_equivalent(theta(mod_2.4), theta(expected))
  ##
  ## Impute HD
  ##
  set.seed(4)
  simdat2 <- augment_rarehomdel(mod_2.4, sb, mb.subsamp, THR)
  if(FALSE){
    saveRDS(simdat2, file=file.path(path, "simdat2.rds"))
  }
  expected <- readRDS(file.path(path, "simdat2.rds"))
  expect_equivalent(simdat2, expected)
  set.seed(5)
  mod_1.4 <- hd4comp(mod_2.4, simdat2, mb.subsamp, mp)
  if(FALSE){
    saveRDS(mod_1.4, file=file.path(path, "final_model.rds"))
  }
  expect_true(mean(rowMaxs(probz(dropSimulated(mod_1.4))) > 0.99) > mean(rowMaxs(probz(sb)) > 0.99))
  expected <- readRDS(file.path(path, "final_model.rds"))
  expect_equivalent(mod_1.4, expected)
  expect_equivalent(theta(mod_1.4), theta(expected))
  expect_equivalent(assays(mod_1.4), assays(expected))
})

test_that("homdeldup_model", {
  library(SummarizedExperiment)
  set.seed(5)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_029"]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  mp <- McmcParams(iter=400, burnin=500)
  mod_1.4 <- homdeldup_model(mb.subsamp, mp)
  expected <- readRDS(file.path(path, "final_model.rds"))
  expect_equivalent(theta(mod_1.4), theta(expected), tolerance=0.05)
})


.test_that("", {
  ##trace(explore_multibatch, browser)
  mb <- revertToMultiBatch(sb)
  restricted <- fit_restricted(mb, sb, THR, model="MBP3",
                               use_restricted=TRUE)
  x <- assays(mb) %>%
    filter(batch==5)
  plot(density(x$oned, adjust=1/2))

  ## potential problems with duplication component
  ## rare duplications can cause problems
  expect_true(any(round(p(restricted), 2)[, 3] < 0.05))
  dat <- filter(assays(sb), oned > 0.1)
  cutoff <- duplication_cutoff(assays(sb))
  restricted2 <- augment_rarecomponent(restricted, sb,
                                      rare_component_restricted=3,
                                      rare_component_sb=4,
                                      diploid_component_sb=3,
                                      use_restricted_theta=TRUE)
  restricted3 <- warmup(assays(restricted2), "MBP3")
  mcmcParams(restricted3) <- mp
  restricted4 <- posteriorSimulation(restricted3)
  full <- mcmcWithHomDel(mb, sb, restricted4, THR)
  gmodel <- genotype_model(full, snpdat)
  expect_identical(mapping(gmodel),
                   as.character(0:3))
})
