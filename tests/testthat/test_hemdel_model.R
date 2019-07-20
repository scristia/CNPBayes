context("hemizygous deletion models")

test_that("summarize CNP_014", {
  set.seed(9209)
  library(panc.data)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_014")
  if(!dir.exists(path)) dir.create(path)
  data(cnp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[14, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
})

test_that("hemdel pipeline", {
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_014")
  set.seed(9209)
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  index <- overlapsAny(snp_se, rowRanges(cnp_se)["CNP_014"])
  snpdat <- snp_se[index, ]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  mp <- McmcParams(iter=400, burnin=500)
  sb <- warmup(assays(mb.subsamp), "SBP2", "SB2")
  mcmcParams(sb) <- mp
  sb <- posteriorSimulation(sb)
  finished <- stop_early(sb)
  expect_false(finished)

  mb <- warmup(assays(mb.subsamp), "MBP2", "MB1")
  expect_error(genotype_model(mb))
  mcmcParams(mb) <- mp
  mb <- posteriorSimulation(mb)
  expect_identical(modelName(mb), "MBP2")
  gmodel <- genotype_model(mb, snpdat)
  path2 <- system.file("extdata", package="cnp.models")
  expected <- readRDS(file.path(path2, "CNP_014.rds"))
  expect_identical(mapping(gmodel), mapping(expected))
  if(FALSE){
    saveRDS(gmodel, file=file.path(path, "CNP_014.rds"))
  }
})

test_that("hemdel_model", {
  set.seed(9209)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_014")
  data(cnp_se, package="panc.data")
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  mp <- McmcParams(iter=400, burnin=500)
  mb <- hemdel_model(mb.subsamp, mp)
  expected <- readRDS(file.path(path, "CNP_014.rds"))
  expect_equivalent(theta(mb), theta(expected))
})

test_that("hemideletion_model", {
  set.seed(9209)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_014")
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  snpdat <- subsetByOverlaps(snp_se, cnp_se[14])
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  summaries(mb.subsamp)$deletion_cutoff <- -0.25
  mp <- McmcParams(iter=400, burnin=500)
  model.list <- hemideletion_models(mb.subsamp, snpdat, mp)
  posthoc <- posthoc_checks(model.list)
  appears_diploid <- not_duplication(model.list[[2]])
  if(appears_diploid){
    model <- model.list[[1]]
  } else {
    ix <- which.min(posthoc$bic)
    model <- model.list[[ix]]
  }
  if(FALSE){
    ggMixture(model.list[[1]]) + xlim(c(-4, 1))
    ggMixture(model.list[[2]]) + xlim(c(-4, 1))
    m2 <- dropSimulated(model.list[[2]])
    ggMixture(m2) + xlim(c(-4, 1))
  }
})
