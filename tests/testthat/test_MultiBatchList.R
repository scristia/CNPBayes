context("MultiBatch classes")

test_that("revised_constructors", {
  library(SummarizedExperiment)
  ##
  ##  - get burnin running for multibatch
  ##
  mb <- MultiBatchModel2()
  mcmcParams(mb) <- McmcParams(iter=iter(mb), thin=thin(mb),
                               burnin=burnin(mb))
  expect_true(validObject(mb))
  mp <- McmcParams(iter=10, burnin=1, nStarts=4)
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mcmcParams(sb) <- mp
  mb2 <- as(sb, "MultiBatch")
  posteriorSimulation(mb2)
  expect_true(validObject(mb2))

  dat <- assays(mb2)
  tmp <- MultiBatch("SB3",
                    data=dat,
                    parameters=parameters(mb2))
  identical(parameters(mb2), parameters(tmp))

  ##
  ## By default, we assume no downsampling
  ##
  expect_identical(assays(mb2), downSampledData(mb2))
  ## down sampling is completely random
  down_sample(mb2) <- sample(seq_len(nrow(mb2)), 1000, replace=TRUE)
  expect_error(validObject(mb2))
  ## safer to do downSampleModel
  mb3 <- downSampleModel(mb2, 1000)
  expect_true(validObject(mb3))
  dataMean(mb3)
  ##
  ## MCMC should be run AFTER downsampling.
  ##
  runBurnin(mb3)
  mb3 <- posteriorSimulation(mb3)
  expect_true(validObject(mb3))
  ddat <- downSampledData(mb2)
  expect_identical(nrow(ddat), 1000L)

  ##
  ## Because down_sample for 'mb2' has been specified, the older class automatically receives the downSampledData
  ##
  sb2 <- as(mb3, "MultiBatchModel")
  expect_true(validObject(sb2))
  expect_identical(z(mb3), z(sb2))
  expect_identical(mu(mb3), mu(sb2))
  expect_identical(theta(chains(mb3)), theta(chains(sb2)))
  ##
  ## Turn one specific model (e.g., MB3) into a list with multiple chains
  sb2.list <- as(mb3, "list")
  expect_is(sb2.list[[1]], "MultiBatchModel")
  ## each element of the list should have nStarts set to 1
  expect_identical(nStarts(sb2.list[[1]]), 1L)

  ##
  ## Convert a list with multiple chains back to one single MultiBatch
  ##
  ## old approach
  ##trace(combine_batch, browser)
  tmp <- combine_batch(sb2.list, batches=1L)
  expect_is(tmp, "MultiBatchModel")
  ch <- chains(tmp)
  ch2 <- combineChains(sb2.list)
  expect_identical(ch, ch2)

  combined.mb <- combineModels(sb2.list)
  expect_false(label_switch(combined.mb))
  ## We've combined the 4 chains into a single chain -- nStarts should be 1
  expect_identical(nStarts(combined.mb), 1L)

  ## accessors for values
  theta(mb2)
  sigma2(mb2)
  sigma(mb2)
  p(mb2)
  nu.0(mb2)
  sigma2.0(mb2)
  mu(mb2)
  tau2(mb2)
  logPrior(mb2)
  log_lik(mb2)
  head(probz(mb2))
  head(u(mb2))
  head(z(mb2))

  ## updates for values
  theta(mb2) <- theta(mb2)
  sigma2(mb2) <- sigma2(mb2)
  sigma(mb2) <- sigma(mb2)
  p(mb2) <- p(mb2)
  nu.0(mb2) <- nu.0(mb2)
  sigma2.0(mb2) <- sigma2.0(mb2)
  mu(mb2) <- mu(mb2)
  tau2(mb2) <- tau2(mb2)
  logPrior(mb2) <- logPrior(mb2)
  log_lik(mb2) <- log_lik(mb2)
  probz(mb2) <- probz(mb2)
  u(mb2) <- u(mb2)
  z(mb2) <- z(mb2)

  chains(mb2) <- chains(mb2)
  parameters(mb2) <- parameters(mb2)
  values(mb2) <- values(mb2)
  flags(mb2)
  flags(mb2) <- flags(mb2)

  ## accessors for parameters
  iter(mb2)
  burnin(mb2)
  thin(mb2)
  iter(mb2) <- 10
  burnin(mb2) <- 10
  thin(mb2) <- 10
  hyperParams(mb2)
  mcmcParams(mb2)
  ## replacement
  hyperParams(mb2) <- hyperParams(mb2)
  mcmcParams(mb2) <- mcmcParams(mb2)
  k(mb2)

  ## flags
  label_switch(mb2)
  label_switch(mb2) <- label_switch(mb2)

  ## summaries
  dataMean(mb2)
  dataPrec(mb2)
  dataMean(mb2) <- dataMean(mb2)
  dataPrec(mb2) <- dataPrec(mb2)
  marginal_lik(mb2)
  marginal_lik(mb2) <- marginal_lik(mb2)

  ## data
  assays(mb2)
  oned(mb2)
  batch(mb2)
  oned(mb2) <- oned(mb2)
  assays(mb2) <- assays(mb2)
})

test_that("starting_values", {
  mb <- MultiBatch()
  expect_true(validObject(mb))

  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  dat <- tibble(oned=y(sb),
                batch=1L) %>%
    mutate(id=seq_len(nrow(.)),
           id=as.character(id))

  mb <- MultiBatch(data=dat)
  mb2 <- MultiBatch(data=dat,
                    down_sample=sort(sample(seq_len(nrow(dat)), 500, replace=TRUE)))
  expect_true(validObject(mb2))
})

test_that("findSurrogates", {
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mb <- as(sb, "MultiBatch")
  assays(mb)[["provisional_batch"]] <- sample(letters[1:15], nrow(mb),
                                              replace=TRUE)
  ##
  ## starting values from modelValues2 are unlikely to be very good
  ##
  mb2 <- findSurrogates(mb)
  iter(mb2) <- 100
  burnin(mb2) <- 50
  thin(mb2) <- 2
  expect_true(validObject(mb2))
  mb2 <- posteriorSimulation(mb2)
})



test_that("downsampling_with_surrogates", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mcmcParams(mb) <- McmcParams(iter=10, thin=2,
                               burnin=50,
                               nStarts=4)
  mb1 <- as(mb, "MultiBatch")
  expect_true(validObject(mb1))
  runBurnin(mb1)
  posteriorSimulation(mb1)

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the chains
  chains(mb1) <- mcmc_chains( specs(mb1), parameters(mb1) )
  mb1 <- posteriorSimulation(mb1)  ## yes

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the summary values
  summaries(mb1) <- modelSummaries( specs(mb1) )
  posteriorSimulation(mb1)  ## yes

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the current values
  summaries(mb1) <- modelSummaries( specs(mb1) )
  posteriorSimulation(mb1)  ## yes

  mp <- McmcParams(iter=10, burnin=10, thin=1, nStarts=4)
  mcmcParams(mb1) <- mp
  ## data mean / precision should have 3 rows!!
  mb2 <- downSampleModel(mb1, 1000)
  mb2 <- posteriorSimulation(mb2)
  expect_true(validObject(mb2))
  ##
  ## simulate plates from the true batches
  ##
  dat <- assays(mb2)
  dat.list <- split(dat, dat$batch)
  plt.list <- list(letters[1:3], letters[4:6], letters[7:9])
  plt.list2 <- map2(dat.list, plt.list,
                    function(x, plt){
                      sample(plt, nrow(x), replace=TRUE)
                    })
  dat$true_batch <- B <- dat$batch
  pb <- rep(NA, nrow(dat))
  pb[ B==1 ] <- plt.list2[[1]]
  pb[ B==2 ] <- plt.list2[[2]]
  pb[ B==3 ] <- plt.list2[[3]]
  dat$provisional_batch <- pb
  dat %>%
    group_by(provisional_batch) %>%
    summarize(n=n(),
              batch=unique(true_batch))
  assays(mb2) <- dat
  expect_true(validObject(mb2))
  mb3 <- findSurrogates(mb2)
  expect_true(validObject(mb3))
  burnin(mb3) <- 0
  mb4 <- posteriorSimulation(mb3)
})

test_that("upsample", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mcmcParams(mb) <- McmcParams(iter=200, thin=2,
                               burnin=200,
                               nStarts=4)
  mb1 <- as(mb, "MultiBatch")
  mb2 <- downSampleModel(mb1, 1000L)
  burnin(mb2) <- 0L
  mb3 <- posteriorSimulation(mb2)
  expect_identical(specs(mb3)$number_obs,
                   specs(mb1)$number_obs)
  expect_identical(names(values(mb3)),
                   names(modes(mb3)))

  ## downsample with specific index
  i <- sort(sample(seq_len(nrow(mb1)), size=1000L, replace=TRUE))
  mb2 <- downSampleModel(mb1, i=i)
  expect_identical(down_sample(mb2), i)
  ##
  ## with upsampling, we would like to infer the probability of the copy number
  ## assignment from the downsampled data.
  ##
  ## Steps:
  ## 1.  evaluate convergence, etc.
  ## 2.  set current values to modal ordinates
  ## 3.  compute z probabilities for full data from modal ordinates
  ##     (all ordinates are modal except for u)
  ## 4.  give z probabilities from step 3, update z in the standard way
  ## Should we automatically upSample?
  mb4 <- upSampleModel(mb3)
  z1 <- z(mb4)
  z2 <- values(mb4)[["z_up"]]
  z2 <- z2[ down_sample(mb4) ]
  ## check that z is at the modal ordinate for the sampled data
  expect_identical(z1, z2)
})

test_that("Pooled model", {
  mb <- MultiBatchP()
  ## with data
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mb1 <- MultiBatchP(data=assays(mb))
  expect_identical(ncol(sigma(mb1)), 1L)

  tmp <- downSampleModel(mb1)
  tmp2 <- as(tmp, "MultiBatchPooled")

  data(MultiBatchPooledExample)
  mbp <- MultiBatchPooledExample
  mp2 <- as(mbp, "MultiBatchP")
  expect_true(validObject(mp2))
  expect_is(mp2, "MultiBatchP")

  iter(mp2) <- 10L
  burnin(mp2) <- 10L
  mp3 <- posteriorSimulation(mp2)
  log_lik(chains(mp3))
})

test_that("Plotting", {
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ch2 <- McmcChains2(mc=chains(mb),
                     iter=iter(mb),
                     k=k(mb),
                     batch=nBatch(mb)) %>%
    as("list")

  ch.list <- listChains(mb)
  expect_identical(ch2, ch.list)
  if(FALSE){
    ggplot(ch2[["theta"]], aes(s, theta, group=b)) +
      geom_line() +
      facet_wrap(~k)
  }
})

test_that("Empty MBL", {
  specs <- modelSpecs()
  data <- modelData()
  down_sample <- seq_len(nrow(data))
  num_models <- nrow(specs)
  iter <- 1000L; thin <- 1L; nStarts <- 4L
  hp <- Hyperparameters()
  mp <- McmcParams(iter=iter,
                   thin=thin,
                   nStarts=nStarts)
  parameters <- modelParameters(hp=hp, mp=mp)
  chains <- listChains2(specs, parameters)
  current_values <- listValues(specs,
                               data[down_sample, , drop=FALSE],
                               hp)
  summaries <- listSummaries(specs)
  flags <- listFlags(specs)
  MultiBatchList()
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
  options(warn=2, error=recover)

  data <- assays(mb1)
  down_sample <- seq_len(nrow(data))
  specs <- modelSpecs(data=data,
                      down_sample=down_sample)
  num_models <- nrow(specs)
  expect_identical(num_models, 8L)
  iter=1000L; thin=1L; nStarts=4L
  hp=Hyperparameters()
  mp=McmcParams(iter=iter, thin=thin, nStarts=nStarts)
  parameters=modelParameters(hp=hp, mp=mp)
  tmp <- listChains2(specs, parameters)
  tmp <- listValues(specs, data[down_sample, , drop=FALSE], hp)
  tmp <- listSummaries(specs)
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
  expect_identical(L, 15L)
  ##
  ## With coercion, we generate a length-1 list
  ##
  mbl <- as(mb1, "MultiBatchList")
  expect_identical(length(mbl), 1L)
  expect_identical(assays(mbl), assays(mb1))
  expect_identical(down_sample(mbl), down_sample(mb1))
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

  expect_identical(seq_along(mbl), 1:15)
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
  ## no code for pooled models yet
  is_pooled <- substr(specs(mbl)$model, 3, 3) == "P"
  mbl <- mbl[ !is_pooled ]
  mbl2 <- posteriorSimulation(mbl)
  expect_true(!any(is.na(theta(chains(mbl2)[[1]]))))
  ##
  ## with downsampling
  ##
  mbl3 <- mbl
  i <- sort(sample(seq_len(nrow(mbl)), replace=TRUE, size=750))
  tmp <- downSampleModel(mbl[[5]], i=i)
  expect_true(validObject(tmp))

  mbl3 <- downSampleModel(mbl, N=750)
  expect_identical(length(u(mbl3[[1]])), 750L)
  mb <- mbl3[[1]]
  expect_true(validObject(mb))
  posteriorSimulation(mb)
  expect_true(all(is.na(theta(chains(mbl3)[[1]]))))
  mbl4 <- posteriorSimulation(mbl3)
  expect_true(all(!is.na(theta(chains(mbl4)[[1]]))))

  mbl5 <- upSampleModel(mbl4)
  pz <- values(mbl5)[[1]][["probz_up"]]
  expect_identical(nrow(pz), specs(mbl5)$number_obs[1])
  expect_false(any(is.na(pz)))
})
