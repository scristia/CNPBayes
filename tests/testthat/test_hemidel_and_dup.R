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
  ##
  ## For hemizygous deletions/duplications, cutoff should depend on batch -- assume it is the same distance from the mode
  ##   - easy to calculate mode of diploid component across batches
  ##   - from single batch model, we have an overall duplication mode that could be shifted by batch
  ##
  library(SummarizedExperiment)
  path <- file.path(system.file("extdata", package="CNPBayes"),
                    "CNP_147")
  data(snp_se, package="panc.data")
  data(cnp_se, package="panc.data")
  g <- rowRanges(cnp_se)["CNP_147"]
  snpdat <- snp_se[overlapsAny(snp_se, g), ]
  mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
  mb <- hemdeldup_model2(mb.subsamp, mp, THR=-0.25)
  mb <- genotype_model(mb, snpdat)
  sb <- genotype_model(sb, snpdat)
  expect_identical(mapping(mb),
                   mapping(sb))
  expect_identical(mapping(mb), c("1", "2", "3"))
  if(FALSE){
    require(grid)
    bafdat <- join_baf_oned(mb, snpdat)
    figs <- list_mixture_plots(mb, bafdat)
    mixture_layout(figs, augmented=TRUE)
  }
})
