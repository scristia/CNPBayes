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

test_that("downsampling with lists", {
  skip("need to fix this unit test")
  ##
  ## with downsampling
  ##
  i <- sort(sample(seq_len(nrow(mbl)), replace=TRUE, size=750))
  ds <- MultiBatchList(data=assays(mb1)[i, ])
  expect_true(validObject(ds))
  iter(ds) <- 100
  burnin(ds) <- 50
  ds <- posteriorSimulation(ds)
  expect_true(validObject(ds))
})

test_that("upsampling with lists", {
  skip("Might be easiest to restrict upsampling to non-list data -- after model selection  ")
  ## Easiest to do upsampling once we have selected a model
  mbl <- full.model
  mbl5 <- upSampleModel(ds, )
  pz <- current_values(mbl5)[[1]][["probz_up"]]
  expect_identical(nrow(pz), specs(mbl5)$number_obs[1])
  expect_false(any(is.na(pz)))
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


test_that("Running multiple chains for multiple models", {
  set.seed(34)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample %>%
    as("MultiBatch") %>%
    "["(sort(sample(seq_len(nrow(.)), 500, replace=TRUE))) %>% {
    MultiBatchList(data=assays(.),
                   parameters=parameters(.))
    }
  orig <- mb[["MB2"]]
  nStarts(mb) <- 1
  iter(mb) <- 0
  burnin(mb) <- 500
  mb2 <- lapply(mb, posteriorSimulation)
  loglik <- sapply(mb2, log_lik) %>%
    sort(decreasing=TRUE) %>%
    head
  mb3 <- mb2[names(loglik)] %>%
    listToMultiBatchList
  ## Might want to specify more iterations here
  nStarts(mb3) <- 3
  iter(mb3) <- 150
  burnin(mb3) <- 300
  for(i in seq_along(mb3)){
    cat(".")
    m <- as(mb3[[i]], "list") %>%
      posteriorSimulation %>%
      combineModels
    fails_gr(m) <- is_high_mpsrf(m)
    mb3[[i]] <- m
  }
  mb3@parameters[["mp"]]@iter <- iter(m)
  mb4 <- mb3[ !sapply(mb3, anyWarnings) ] %>%
    compute_marginal_lik
  ## For marginal likelihood
  x <- mb4[ order(marginal_lik(mb4), decreasing=TRUE) ]
  ##
  ## reset the initial model
  ##
  mb2 <- reset(orig, x[[1]])
  expect_true(validObject(mb2))
  expect_identical( specs(mb2), specs(x[[1]]))
  ##expect_equivalent(summaries2(mb2), summaries2(x[[1]]))
  expect_identical(current_values2(mb2), current_values2(x[[1]]))
})


test_that("augment data for MultiBatch", {
  set.seed(123)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mb <- posteriorSimulation(mb)
  if(FALSE){
    ch.list <- as(chains(mb), "list")
    ggplot(ch.list[["theta"]], aes(s, value, group=b)) +
      geom_line(aes(color=b)) +
      facet_wrap(~k)
  }
  filter <- dplyr::filter
  ## remove homozygous deletions from batch 3
  ## - not a very realistic example because this is not a rare deletion
  dat <- assays(mb) %>%
    filter(! (batch == 3 & oned < -0.5) )
  ## we have problems with this because one batch has only 2 components
  mb1 <- MultiBatch(data=dat)
  iter(mb1) <- 250L
  burnin(mb1) <- 100L
  mb1 <- posteriorSimulation(mb1)
  if(FALSE){
    ## first and second component completely overlap for batch 3
    ggMixture(mb1)
  }
  mbl <- MultiBatchList(data=assays(mb1),
                        mp=mcmcParams(mb1))
  ##
  ## - check for batches with few observations
  ## - augment data if needed
  ## - fit as per usual
  ##
  mbl2 <- augmentData2(mbl)
  expect_true(validObject(mbl2))
  ## 10 observations added
  expect_equal(nrow(mbl2), 1451L, tolerance=3)
  mb3 <- mbl2[[ which(specs(mbl2)$model == "MB3") ]]
  mb3 <- posteriorSimulation(mb3)
  if(FALSE) ggMixture(mb3)
})

test_that("fix probz mcmc2", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  iter(mb) <- 20
  burnin(mb) <- 20
  nStarts(mb) <- 2
  max_burnin(mb) <- 20
  mb2 <- mcmc2(mb)
  expect_true(all(rowSums(probz(mb2))==1))
})

##
## This test is slow
##
test_that("singleBatchGuided", {
  set.seed(1968)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb), burnin=100L,
                        iter=150L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 100L
  sb <- mcmc2(mbl[[ "SB3" ]])
  ##expect_true( convergence(sb) )

  x <- mbl[ modelName(mbl) != "SB3" ]
  xx <- singleBatchGuided(x, sb)
  expect_identical(names(xx), modelName(x))
  expect_is(xx, "list")
  expect_is(xx[[1]], "list")
  expect_identical(length(xx[[1]]), nStarts(x))
  expect_identical(length(xx), length(x))
  x <- mbl[["MB3"]]
  xx <- singleBatchGuided(x, sb)
  th <- theta(sb)
  expect_equal(colMeans(theta(xx)), th[1, ], tolerance=0.15)
})

test_that("what is wrong with pooled-variance models", {
  set.seed(2000)
  truth <- simulateData(N = 1000,
                        theta = c(-2, -0.4, 0),
                        sds = c(0.05, 0.05, 0.05),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10),
                        df=100)
  mp <- McmcParams(iter = 100, burnin = 100)
  hp <- Hyperparameters(k=3)
  m <- SingleBatchPooled(y(truth), hp, mp)
  model <- posteriorSimulation(m)
  expect_equal(sigma(model)[[1]], 0.05, tolerance=0.02)
  expect_equal(theta(model)[1, ], theta(truth)[1, ], tolerance=0.05)
  if(FALSE) ggMixture(model)

  mb <- as(m, "MultiBatchP")
  mb2 <- posteriorSimulation(mb)
  expect_equal(sigma(mb2)[, 1], 0.05, tolerance=0.02)
  expect_equal(theta(mb2)[1, ], theta(truth)[1, ], tolerance=0.1)
  if(FALSE) ggMixture(mb2)

  ## pooled model with 3 batches
  set.seed(1438)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  sbp <- as(MultiBatchModelExample, "MultiBatchPooled")
  expect_identical(length(sigma(sbp)), 3L)

  mb <- as(MultiBatchModelExample, "MultiBatchP")
  expect_identical(dim(sigma(mb)), c(3L, 1L))
  nStarts(mb) <- 1L
  mb2 <- posteriorSimulation(mb)
  if(FALSE) ggMixture(mb2)

  ## Try guided function
  set.seed(150)
  mbl <- MultiBatchList(data=assays(mb), burnin=200L,
                        iter=300L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 200L
  sb1 <- mbl[[ "SB3" ]]
  sb <- posteriorSimulation(sb1)
  ml <- mcmcList(sb)
  effsize <- ml %>%
    effectiveSize %>%
    median
  gelman_rub <- gelman_rubin(ml)$mpsrf
  flags(sb)[["small_effsize"]] <- effsize < min_effsize(sb)
  flags(sb)[["fails_GR"]] <- gelman_rub > min_GR(sb)
  smallcount <- flags(sb)$.internal.counter > 10
  fl <- c(label_switch=flags(sb)[["label_switch"]],
          small_effsize=flags(sb)[["small_effsize"]],
          fails_GR=flags(sb)[["fails_GR"]],
          high_internal_counter=smallcount)
  any_flags <- any(fl)
  expect_identical( convergence(sb), !any_flags )
  x <- mbl[ modelName(mbl) != "SB3" ]
  ##
  ## Something is wrong with the pooled-variance models
  ## With guided approach, we first fit SB3 (above).
  ## Since SB3 converged, we then fit SBP3, MB3, and MBP3 models
  ## with starting values near the SB3 parameters
  ##
  xx <- singleBatchGuided(x, sb)
  sbp <- xx[["SBP3"]][[1]]
  iter(sbp) <- 250L
  burnin(sbp) <- 100L
  sbp <- posteriorSimulation(sbp)
  th <- theta(sbp) %>%
    as.numeric
  expect_equal(th, theta(sb)[1, ], tolerance=0.05)
  if(FALSE) ggMixture(sbp)
  ##
  ## MultiBatch
  ##
  mb <- xx[["MB3"]][[1]]
  mcmcParams(mb) <- mcmcParams(sbp)
  ## This should be faster because we are already in region of high posterior probability
  mb <- posteriorSimulation(mb)
  if(FALSE) ggMixture(mb)
  expect_equal(colMeans(theta(mb)), theta(sb)[1, ],
               tolerance=0.05)
  ##
  ## Guided MCMC2
  ## Now we try the guided version of mcmc2, which basically implements the above steps
  ##mbl <- mbl[ names(mbl) != "SB3" ]
  ##mbl2 <- mcmc2(mbl, sb)
})

test_that("MultiBatchP <-> MultiBatchPooled", {
  x <- MultiBatchPooledExample
  x2 <- as(x, "MultiBatchP") %>%
    as("MultiBatchPooled")
  expect_true(all(names(modes(x)) %in% names(modes(x2))))
})

test_that("augment data from pooled model", {
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb), burnin=100L,
                        iter=300L, max_burnin=100L)
  mbp <- mbl[specs(mbl)$model == "MBP3"]
  mbp2 <- toSingleBatch(mbp[[1]])
  expect_true(validObject(mbp2))
  mb <- mbl[specs(mbl)$model == "MB3"]
  mb2 <- toSingleBatch(mb[[1]])
  expect_true(validObject(mb2))
  ##
  ## We need to be able to create a model with a single batch in order to assess whether there are homozygous deletions in the marginal distribution
  ##
  ##  - to do this, various aspects of the model must be reparameterized to have the appropriate dimensions
  ##
  mbp2 <- augmentData2(mbp)
})

test_that("Data not in batch-order", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ix <- sample(seq_len(nrow(mb)), 500, replace=FALSE)
  ##
  ## object not instantiated because data is not in batch order
  ##
  expect_error(MultiBatch(data=assays(mb)[ix, ]))
  expect_error(MultiBatchP(data=assays(mb)[ix, ]))
  expect_error(MultiBatchList(data=assays(mb)[ix, ]))
  mb <- MultiBatchModelExample
  ##
  ## No error is generated with the old classes
  ##
  batch(mb) <- sample(batch(mb), 1500, replace=FALSE)
  validObject(mb)
  expect_error(as(mb, "MultiBatch"))
  expect_error(as(mb, "MultiBatchP"))
  expect_error(as(mb, "MultiBatchList"))
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

