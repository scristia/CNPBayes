context("CNV model selection")
## CNP_001

test_that("cnv pipeline", {
  library(SummarizedExperiment)
  path <- system.file("extdata", package="CNPBayes")
  cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_001", ]
  snp_se <- readRDS(file.path(path1, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se)
  mb <- readRDS(file.path(path, "mb_subsamp.rds"))
  model <- cnv_models(mb, rowRanges(cnp_se), snp_se)
  expect_identical(modelName(model), "MBP3")
  expect_identical(mapping(model), c("0", "2", "2"))
  if(FALSE){
    ## verify probabilities well calibrated with constrained mcmc
    maxpz <- tibble(prob=rowMaxs(probz(model)),
                    batch=batch(model),
                    z=map_z(model))
    cutoffs <- group_by(maxpz, z) %>%
      summarize(upper_quartile=quantile(prob, 0.75))
    ggplot(maxpz, aes(z, prob)) +
      geom_jitter(width=0.1) +
      facet_wrap(~batch)
    bafdat <- join_baf_oned(model, snp_se)
    figs <- list_mixture_plots(model, bafdat)
    mixture_layout(figs, augmented=TRUE)
  }
})

##
## Common deletion
##
test_that("common deletion", {
  path <- system.file("extdata", package="CNPBayes")
  set.seed(555)
  cnp_se <- readRDS(file.path(path, "cnp_se.rds"))["CNP_022", ]
  snp_se <- readRDS(file.path(path, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se)
  path2 <- file.path(path, "CNP_022")
  mb.subsamp <- readRDS(file.path(path2, "mb_subsamp.rds"))
  mp <- McmcParams(iter=400, burnin=500)
  model.list <- deletion_models(mb.subsamp, snp_se, mp)
  model <- choose_model(model.list, mb.subsamp)
  expect_identical(mapping(model), c("0", "1", "2"))
})

## homozygous deletion and duplication
test_that("homdeldup_model", {
  library(SummarizedExperiment)
  path1 <- system.file("extdata", package="CNPBayes")
  path <- file.path(path1, "CNP_029")
  cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_029", ]
  snp_se <- readRDS(file.path(path1, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se)
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  mp <- McmcParams(iter=400, burnin=500)
  set.seed(5)
  model <- homdeldup_model(mb.subsamp, mp)
  gmodel <- genotype_model(model, snp_se)
  expect_identical(mapping(gmodel), c("0", "1", "2", "3"))
  if(FALSE){
    model2 <- select_highconfidence_samples(model, snp_se)
    bafdat <- join_baf_oned(model, snp_se)
    figs <- list_mixture_plots(model, bafdat)
    mixture_layout(figs, augmented=TRUE)    
  }
})

##
## hemizygous deletion
##
test_that("hemideletion_model", {
  set.seed(9209)
  path1 <- system.file("extdata", package="CNPBayes")
  path <- file.path(path1, "CNP_014")
  cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))
  snp_se <- readRDS(file.path(path1, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se["CNP_014"])
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  summaries(mb.subsamp)$deletion_cutoff <- -0.25
  mp <- McmcParams(iter=400, burnin=500)
  model.list <- hemideletion_models(mb.subsamp, snp_se, mp)
  model <- choose_model(model.list)
  expect_identical(modelName(model), "MBP2")
})

test_that("Hemizygous deletion +/- duplication pipeline", {
  ##
  ## For hemizygous deletions/duplications, cutoff should depend on batch -- assume it is the same distance from the mode
  ##   - easy to calculate mode of diploid component across batches
  ##   - from single batch model, we have an overall duplication mode that could be shifted by batch
  ##
  set.seed(9209)
  path1 <- system.file("extdata", package="CNPBayes")
  path <- file.path(path1, "CNP_147")
  cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_147", ]
  snp_se <- readRDS(file.path(path1, "snp_se.rds"))
  snp_se <- subsetByOverlaps(snp_se, cnp_se)
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  summaries(mb.subsamp)$deletion_cutoff <- -0.25
  mp <- McmcParams(iter=400, burnin=500)
  set.seed(9209)
  model.list <- hemideletion_models(mb.subsamp, snp_se, mp)
  model <- choose_model(model.list)
  expect_identical(mapping(model), c("2", "3"))
})


