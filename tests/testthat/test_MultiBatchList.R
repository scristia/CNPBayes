contect("MultiBatch classes")

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
  expect_false(validObject(mb2))
  ## safer to do downSampleModel
  mb3 <- downSampleModel(mb2)
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
  expect_identical(nStarts(combined.mb), length(sb2.list))

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
  label_switch(mb2) <- label_switch(mb2)a

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
  iter(mb2) <- 200
  burnin(mb2) <- 100
  expect_true(validObject(mb2))
  mb2 <- posteriorSimulation(mb2)
})



test_that("downsampling", {
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
  mb2 <- downSampleModel(mb1)
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
