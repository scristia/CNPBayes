context("MultiBatchList classes")

test_that("Empty MBL", {
  specs <- modelSpecs()
  data <- modelData()
  num_models <- nrow(specs)
  iter <- 1000L; thin <- 1L; nStarts <- 4L
  hp <- Hyperparameters()
  mp <- McmcParams(iter=iter,
                   thin=thin,
                   nStarts=nStarts)
  parameters <- modelParameters(hp=hp, mp=mp)
  chains <- listChains2(specs, parameters)
  current_values <- listValues(specs,
                               data,
                               hp)
  summaries <- listSummaries(specs)
  flags <- listFlags(specs)
  expect_true(validObject(MultiBatchList()))
})

test_that("MBL with data", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mcmcParams(sb) <- McmcParams(iter=200, thin=2,
                               burnin=200,
                               nStarts=4)
  mb1 <- as(sb, "MultiBatch")
  ## When a single model is not specified,
  ## generate all available models
  mbl <- MultiBatchList(data=assays(mb1))
  expect_true(validObject(mbl))

  ## data with multiple batches
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb1 <- as(mb, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb1))
  L <- length(mbl)
  expect_identical(L, 14L)
  ##
  ## With coercion, we generate a length-1 list
  ##
  mbl <- as(mb1, "MultiBatchList")
  expect_identical(length(mbl), 1L)
  expect_identical(assays(mbl), assays(mb1))
  ##
  ## MultiBatchList -> MultiBatch
  ##  -- only takes first element from list
  mb <- as(mbl, "MultiBatchList")
  expect_true(validObject(mb))

  mbl <- MultiBatchList(data=assays(mb1))
  ## grab a subset of the models
  mbl2 <- mbl[c(1 , 3, 5)]
  expect_identical(specs(mbl2)$model,
                   specs(mbl)$model[c(1, 3, 5)])

  expect_identical(seq_along(mbl), 1:14)
})

test_that("posterior simulation with MBL", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb1 <- as(mb, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb1))
  ## this involves resizing the chains
  mcmcParams(mbl) <- McmcParams(iter=10L,
                                thin=1L,
                                burnin=5L,
                                nStarts=3L)
  expect_true(validObject(mbl))
  ## no code for pooled models yet
  is_pooled <- substr(specs(mbl)$model, 3, 3) == "P"
  mbl <- mbl[ !is_pooled ]
  mbl2 <- posteriorSimulation(mbl)
  expect_true(!any(is.na(theta(chains(mbl2)[[1]]))))
})

test_that("Running multiple chains for one model", {
  set.seed(34)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample %>%
    as("MultiBatch")
  nStarts(mb) <- 3
  iter(mb) <- 150
  burnin(mb) <- 400
  mbl <- as(mb, "list")
  lapply(mbl, theta)
  lapply(mbl, sigma)
  lapply(mbl, p)
  mbl <- lapply(mbl, posteriorSimulation)
  mb <- combineModels(mbl)
  gr <- mcmcList(mbl) %>%
    gelman_rubin
  gr$mpsrf
  flags(mb)$fails_GR <- gr$mpsrf > 1.2
  expect_true(!anyWarnings(mb))
})



test_that("MultiBatchP <-> MultiBatchPooled", {
  x <- MultiBatchPooledExample
  x2 <- as(x, "MultiBatchP") %>%
    as("MultiBatchPooled")
  expect_true(all(names(modes(x)) %in% names(modes(x2))))
})

test_that("Generate a list of candidate copy number models", {
  mb <- MultiBatchModelExample
  mb <- as(mb, "MultiBatch")
  mapping(mb) <- seq_len(k(mb))
  expect_identical(mapping(mb), seq_len(k(mb)))
  clist <- CnList(mb)
  expect_is(clist, "CnList")
  mb <- clist[[1]]
  expect_is(mb, "MultiBatch")
  expect_identical(mapping(mb), c("0", "1", "2"))
  expect_identical(cmap(mb), "012")
})

