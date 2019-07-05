context("MultiBatch classes")

.test_that <- function(nm, expr) NULL

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

  ##sb <- updateObject(sb)
  mcmcParams(sb) <- mp
  expect_true(validObject(sb))
  mb2 <- as(sb, "MultiBatch")
  mb2 <- posteriorSimulation(mb2)

  iter(mb2) <- 0
  burnin(mb2) <- 25
  mb2 <- posteriorSimulation(mb2)
  nStarts(mb2) <- 3

  mb3 <- as(mb2, "list")
  posteriorSimulation(mb3[[1]])

  mbp <- as(mb2, "MultiBatchP")
  mb.list <- as(mbp, "list")
  lapply(mb.list, posteriorSimulation)
  iter(mbp) <- 10
  mb.list <- as(mbp, "list") %>%
    lapply(posteriorSimulation)

  expect_true(validObject(mb2))
  m <- computeModes(mb2)
  dat <- assays(mb2)
  tmp <- MultiBatch("SB3",
                    data=dat,
                    parameters=parameters(mb2))
  identical(parameters(mb2), parameters(tmp))
})

test_that("list of models with independent starting values", {
  skip("replicateMultiBatch not implemented")
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mb <- sb %>%
    as("MultiBatch")
  mb.list <- replicateMultiBatch(mb)

  replicateMultiBatch <- function(object){


  }
  mb <- MultiBatchModel2()
  mcmcParams(mb) <- McmcParams(iter=iter(mb), thin=thin(mb),
                               burnin=burnin(mb))
})

test_that("no downsampling", {
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mb <- sb %>%
    as("MultiBatch")
  mcmcParams(mb) <- McmcParams(iter=10L, thin=0L,
                               burnin=0L)
  mb2 <- posteriorSimulation(mb)
  ##
  ## By default, we assume no downsampling
  ##
  expect_true(validObject(mb2))
  dataMean(mb2)
  ##
  ## MCMC should be run AFTER downsampling.
  ##
  ix <- sample(seq_len(nrow(mb2)), 1000, replace=TRUE)
  mb3 <- mb2[ix, ]
  runBurnin(mb3)

  mb3 <- posteriorSimulation(mb3)
  expect_true(validObject(mb3))
  ddat <- assays(mb3)
  expect_identical(nrow(ddat), 1000L)

  ##
  ## Because down_sample for 'mb2' has been specified, the older class automatically receives the downSampledData
  ##
  sb2 <- as(mb3, "MultiBatchModel")
  expect_true(validObject(sb2))
  expect_identical(z(mb3), z(sb2))
  expect_identical(mu(mb3), mu(sb2))
  expect_identical(theta(chains(mb3)), theta(chains(sb2)))
})

test_that("working with lists of models", {
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mb <- sb %>%
    as("MultiBatch")
  mcmcParams(mb) <- McmcParams(iter=10L, thin=0L,
                               burnin=0L,
                               nStarts=4)
  ##
  ## Turn one specific model (e.g., MB3) into a list with multiple chains
  mb.list <- as(mb, "list")
  expect_is(mb.list[[1]], "MultiBatch")
  ## each element of the list should have nStarts set to 1
  expect_identical(nStarts(mb.list[[1]]), 1L)
  expect_identical(length(mb.list), nStarts(mb))
  ##
  ## Each model should have a different start simulated from hyperpriors
  ##
  ##  is_true <- all.equal(theta(mb.list[[1]]), theta(mb.list[[2]])) %>%
  ##    isTRUE
  ##  expect_false(is_true)
  ##
  ## Convert a list with multiple chains back to one single MultiBatch
  ##
  ## old approach
  ##trace(combine_batch, browser)
  ##tmp <- combine_batch(mb.list, batches=1L)
  tmp <- combineModels(mb.list)
  expect_is(tmp, "MultiBatch")
  ch <- chains(tmp)
  ch2 <- combineChains(mb.list)
  expect_identical(ch, ch2)

  combined.mb <- combineModels(mb.list)
  expect_false(label_switch(combined.mb))
  ## We've combined the 4 chains into a single chain -- nStarts should be 1
  expect_identical(nStarts(combined.mb), 1L)
})

test_that("starting_values", {
  mb <- MultiBatchP()
  mb <- MultiBatch()
  expect_true(validObject(mb))

  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  dat <- tibble(oned=y(sb),
                batch=1L) %>%
    mutate(id=seq_len(nrow(.)),
           id=as.character(id))
})

test_that("findSurrogates", {
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  ##sb <- updateObject(SingleBatchModelExample)
  mb <- as(sb, "MultiBatch")
  assays(mb)[["provisional_batch"]] <- sample(letters[1:15], nrow(mb),
                                              replace=TRUE)
  ##
  ## starting values from modelValues2 are unlikely to be very good
  ##
  mb2 <- findSurrogates(mb)
  expect_true(validObject(mb2))
  iter(mb2) <- 100
  burnin(mb2) <- 50
  thin(mb2) <- 2
  expect_true(validObject(mb2))
  mb2 <- posteriorSimulation(mb2)
  expect_true(validObject(mb2))
})

test_that("downsampling_with_surrogates", {
  set.seed(456)
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
  ix <- sort(sample(seq_len(nrow(mb1)), 1000, replace=TRUE))
  mb2 <- mb1[ix]
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
  expect_true(validObject(mb4))
})

test_that("upsample", {
  skip('skip upsample tests')
  library(SummarizedExperiment)
  library(tidyverse)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mcmcParams(mb) <- McmcParams(iter=200, thin=2,
                               burnin=200,
                               nStarts=4)
  mb1 <- as(mb, "MultiBatch")
  ix <- sort(sample(seq_len(nrow(mb1)), 1000L, replace=TRUE))
  downsampled.model <- mb1[ix]
  burnin(downsampled.model) <- 0L
  downsampled.model <- posteriorSimulation(downsampled.model)
  expect_identical(names(current_values(downsampled.model)),
                   names(modes(downsampled.model)))

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
  full.model <- mb1
  mb4 <- upSampleModel(downsampled.model, full.model)
  if(FALSE) ggMixture(mb4)
})

test_that("Pooled model", {
  mb <- MultiBatchP()
  ## with data
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb <- as(mb, "MultiBatch")
  mb1 <- MultiBatchP(data=assays(mb))
  expect_identical(ncol(sigma_(mb1)), 1L)
  expect_identical(ncol(sigma_(chains(mb1))), 3L)

  data(MultiBatchPooledExample)
  mbp <- MultiBatchPooledExample
  expect_identical(ncol(sigma2(chains(mbp))), numBatch(mbp))
  mp2 <- as(mbp, "MultiBatchP")
  expect_identical(ncol(sigma_(chains(mp2))), 3L)
  expect_identical(nrow(sigma_(mp2)), 3L)
  expect_true(validObject(mp2))
  expect_is(mp2, "MultiBatchP")

  iter(mp2) <- 100L
  expect_identical(ncol(sigma_(chains(mp2))), 3L)
  burnin(mp2) <- 100L
  expect_identical(ncol(sigma_(chains(mp2))), 3L)
  mp3 <- posteriorSimulation(mp2)
  mp <- as(mp3, "MultiBatchPooled")
  ll <- loglik_multibatch_pvar(mp)
  expect_true(!is.na(ll))

  log_lik(chains(mp3))
  expect_true(validObject(mp3))

  set.seed(321)
  ##
  ## These models should fit the data well because
  ## we are starting from region with high posterior probability
  mb1 <- as(mb, "MultiBatchP")
  mbp <- as(mb1, "MultiBatchPooled")
  tmp <- posteriorSimulation(mbp)
  expect_equal(theta(tmp), theta(mb), tolerance=0.02)
  mb1 <- as(mbp, "MultiBatchP")
  tmp <- posteriorSimulation(mb1)
  expect_equal(theta(tmp), theta(mb), tolerance=0.02)
  ##
  ## Random starts
  ##
  mb1 <- MultiBatch(data=assays(mb1),
                    iter=300,
                    burnin=200)
  ## moves to region of high posterior probability quickly
  tmp <- posteriorSimulation(mb1)
  expect_equal(theta(tmp), theta(mb), tolerance=0.02)
  ## pooled model does not move to region of high posterior prob.
  mb1 <- MultiBatchP(data=assays(mb1),
                     iter=300,
                     burnin=200)
  ## does NOT move to region of high posterior probability quickly
  tmp <- posteriorSimulation(mb1)
  ##expect_true(all.equal(theta(tmp), theta(mb), tolerance=0.02))
})

test_that("Plotting", {
  library(tidyverse)
  data(MultiBatchModelExample)
  mbe <- MultiBatchModelExample
  mb <- as(mbe, "MultiBatch")
  ch <- chains(mb)
  ch.list <- as(ch, "list")
  if(FALSE){
    ggplot(ch.list[["theta"]], aes(s, value, group=b)) +
      geom_line(aes(color=b)) +
      facet_wrap(~k)
  }
})

test_that("predictive", {
  library(tidyverse)
  library(magrittr)
  data(MultiBatchModelExample)
  mbe <- MultiBatchModelExample
  prob <- p(mbe)
  z <- seq_len(k(mbe))
  set.seed(123)
  sample(z, replace=TRUE, size=3, prob=prob[1, ])
  set.seed(123)
  tmp=replicate(1000, sample_components(z, 3, prob[1, ])) %>%
    as.numeric
  phat <- table(tmp)/length(tmp)
  phat <- as.numeric(phat)
  expect_equal(phat, prob[1, ], tolerance=0.01)

  set.seed(123)
  update_predictive(mbe)
  ## simulate from rst
  uu <- u(mbe)
  x <- rst(1)

  set.seed(123)
  tmp <- rlocScale_t(1000, 0, 1, 100, uu[1])
  set.seed(123)
  tmp2 <- rst(1000, uu[1], 100, 0, 1)
  expect_identical(tmp, tmp2)

  ch <- chains(mbe)
  predictive(ch) <- matrix(as.numeric(NA),
                           nrow=iter(mbe),
                           ncol=k(mbe)*nBatch(mbe))
  ##
  ## infrastructure set up, but values not quite right
  ##
  ## shouldn't yhat be the same dimension as y?
  ##  i.e., the probability a value is simulated from batch y depends on the representation of the batches.
  chains(mbe) <- ch
  set.seed(123)
  tmp <- cpp_mcmc(mbe)
  ystar <- predictive(tmp)
  ##ystar2 <- posteriorPredictive(mbe)
  ##qqplot(as.numeric(ystar), ystar2$y)

  mbe2 <- posteriorSimulation(mbe)
  pmeans <- colMeans(p(chains(mbe2)))
  ystar <- predictive(mbe2)

  mbp <- MultiBatchPooledExample
  expect_true(validObject(mbp))
  expect_identical(ncol(sigma_(chains(mbp))), nBatch(mbp))
  expect_identical(ncol(theta(chains(mbp))), nBatch(mbp)*k(mbp))

  tmp <- update_predictiveP(mbp)
  mbp3 <- posteriorSimulation(mbp)
  pred <- predictive(mbp3)
  if(FALSE){
    hist(pred, breaks=500)
    hist(y(mbp3), breaks=200)
  }
  p1 <- mean(y(mbp3) < -1)
  p2 <- mean(pred < -1)
  expect_equal(p1, p2, tolerance=0.05)
})

test_that("likelihood for pooled variance model", {
  mbp <- MultiBatchPooledExample
  logl <- loglik_multibatch_pvar(mbp)
  expect_equal(logl, -124, tolerance=1)
})


