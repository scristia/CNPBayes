context("Hemizygous deletion and duplication")

test_that("summarize CNP_147", {
  library(panc.data)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_147")
  if(!dir.exists(path)) dir.create(path)
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  plates <- colData(cnp_se)$Sample.Plate
  mb.subsamp <- summarize_region(cnp_se[147, ],
                                 provisional_batch=plates,
                                 THR=-1)
  if(FALSE){
    saveRDS(mb.subsamp, file=file.path(path, "mb_subsamp.rds"))
  }
  expected <- readRDS(file.path(path, "mb_subsamp.rds"))
  expect_equivalent(mb.subsamp, expected)
})

test_that("Hemizygous deletion / duplication pipeline", {
  library(SummarizedExperiment)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_029")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_029"]
  snpdat <- snp_se[overlapsAny(snp_se, g), ]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
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
  mb <- explore_multibatch(sb, simdat, THR)
  mb <- genotype_model(mb, snpdat)
  sb <- genotype_model(sb, snpdat)
  expect_identical(mapping(mb),
                   mapping(sb))
  expect_identical(mapping(mb), c("1", "2", "3"))
  gmodel <- mb
  if(FALSE){
    require(grid)
    bafdat <- join_baf_oned(mb, snpdat)
    figs <- list_mixture_plots(mb, bafdat)
    mixture_layout(figs, augmented=TRUE)
  }
})
