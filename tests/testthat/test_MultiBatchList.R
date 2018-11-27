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
  expect_true(validObject(mb2))
  m <- computeModes(mb2)
  dat <- assays(mb2)
  tmp <- MultiBatch("SB3",
                    data=dat,
                    parameters=parameters(mb2))
  identical(parameters(mb2), parameters(tmp))
})

test_that("no downsampling", {
  mb <- SingleBatchModelExample %>%
    as("MultiBatch")
  mcmcParams(mb) <- McmcParams(iter=10L, thin=0L,
                               burnin=0L)
  mb2 <- posteriorSimulation(mb)
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
})

test_that("working with lists", {
  mb <- SingleBatchModelExample %>%
    as("MultiBatch")
  mcmcParams(mb) <- McmcParams(iter=10L, thin=0L,
                               burnin=0L,
                               nStarts=4)
  ##
  ## Turn one specific model (e.g., MB3) into a list with multiple chains
  mb.list <- as(mb, "list")
  expect_is(mb.list[[1]], "MultiBatchModel")
  ## each element of the list should have nStarts set to 1
  expect_identical(nStarts(mb.list[[1]]), 1L)
  expect_identical(length(mb.list), nStarts(mb))
  ##
  ## Each model should have a different start simulated from hyperpriors
  ##
  is_true <- all.equal(theta(mb.list[[1]]), theta(mb.list[[2]])) %>%
    isTRUE
  expect_false(is_true)
  ##
  ## Convert a list with multiple chains back to one single MultiBatch
  ##
  ## old approach
  ##trace(combine_batch, browser)
  tmp <- combine_batch(mb.list, batches=1L)
  expect_is(tmp, "MultiBatchModel")
  ch <- chains(tmp)
  ch2 <- combineChains(mb.list)
  expect_identical(ch, ch2)

  combined.mb <- combineModels(mb.list)
  expect_false(label_switch(combined.mb))
  ## We've combined the 4 chains into a single chain -- nStarts should be 1
  expect_identical(nStarts(combined.mb), 1L)
})

test_that("accessors", {
  mb <- SingleBatchModelExample %>%
    as("MultiBatch")
  mcmcParams(mb) <- McmcParams(iter=10L, thin=0L,
                               burnin=0L,
                               nStarts=4)
  mb2 <- posteriorSimulation(mb)
  ## accessors for values
  theta(mb2)
  sigma2(mb2)
  sigma_(mb2)
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
  sigma_(mb2) <- sigma_(mb2)
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
  current_values(mb2) <- current_values(mb2)
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
  ##assays(mb2) <- assays(mb2)
  ##expect_true(validObject(mb2))
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
                    down_sample=sort(sample(seq_len(nrow(dat)),
                                            500, replace=TRUE)))
  expect_true(validObject(mb2))
})

test_that("findSurrogates", {
  data(SingleBatchModelExample)
  sb <- updateObject(SingleBatchModelExample)
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
  expect_identical(names(current_values(mb3)),
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
  z2 <- current_values(mb4)[["z_up"]]
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
  expect_identical(ncol(sigma_(mb1)), 1L)

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
  ##expect_false(all.equal(theta(tmp), theta(mb), tolerance=0.02))
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
  sample(z, replace=TRUE, size=3, prob=prob)
  set.seed(123)
  tmp=replicate(1000, sample_components(z, 3, prob)) %>%
    as.numeric
  phat <- table(tmp)/length(tmp)
  phat <- as.numeric(phat)
  expect_equal(phat, prob, tolerance=0.01)

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
  pz <- current_values(mbl5)[[1]][["probz_up"]]
  expect_identical(nrow(pz), specs(mbl5)$number_obs[1])
  expect_false(any(is.na(pz)))
})

## this is a much simpler replacement for gibbs
test_that("mcmc2 for MultiBatch", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  set.seed(1234)
  mb2 <- as(mb, "MultiBatch")
  iter(mb2) <- 100
  burnin(mb2) <- 100
  max_burnin(mb2) <- 200
  nStarts(mb2) <- 4
  ##trace(startingValues2, browser)
  mb.list <- startingValues2(mb2)
  mb3 <- mcmc2(mb2)
  ## select tolerance of 3
  ## - ml should be more stable with more iterations
  expect_equal(-341, marginal_lik(mb3)[[1]], tolerance=3)
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
  expect_identical(seq_len(nrow(mb3)), down_sample(mb3))
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
  mbl <- MultiBatchList(data=assays(mb), burnin=200L,
                        iter=300L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 200L
  sb <- mcmc2(mbl[[ "SB3" ]])
  expect_true( convergence(sb) )

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
  mbl <- MultiBatchList(data=assays(mb), burnin=200L,
                        iter=300L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 200L
  sb <- mcmc2(mbl[[ "SB3" ]])
  expect_true( convergence(sb) )
  x <- mbl[ modelName(mbl) != "SB3" ]
  ##
  ## Something is wrong with the pooled-variance models
  ## With guided approach, we first fit SB3 (above).
  ## Since SB3 converged, we then fit SBP3, MB3, and MBP3 models
  ## with starting values near the SB3 parameters
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

.test_that <- function(nm, expr) NULL

test_that("Smarter MCMC for MultiBatchList", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ##
  ## beginning of mcmc2
  ##
  object <- MultiBatchList(data=assays(mb), burnin=100L,
                           iter=300L, max_burnin=100L)
  expect_identical(length(unique(batch(object[["SB3"]]))), 1L)
  expect_identical(length(unique(batch(object[["MB3"]]))), 3L)
  sb <- object[["SBP3"]]
  s <- sigma(sb)
  expect_identical(dim(s), c(1L, 1L))
  sb <- object[["MBP3"]]
  expect_is(sb, "MultiBatchP")
  s <- sigma(sb)
  ch <- sigma(chains(sb))
  expect_identical(dim(s), c(3L, 1L))
  expect_identical(dim(summaries(sb)[["data.prec"]]), c(3L, 1L))
  expect_identical(dim(ch), c(iter(sb), nBatch(sb)))

  x <- object[["SBP3"]]
  expect_is(x, "MultiBatchP")
  expect_true(validObject(x))
  x <- object[["MBP3"]]
  expect_is(x, "MultiBatchP")
  expect_true(validObject(x))

  mlist <- listModelsByDecreasingK(object)
  expect_identical(names(mlist), as.character(4:1))
  s <- sigma(mlist[["3"]][["MBP3"]])
  expect_identical(dim(s), c(3L, 1L))
  ##k4 <- fitModelK(mlist[[1]])
  k3 <- fitModelK(mlist[[2]])
  expect_is(k3, "MultiBatchList")
  ## models in list ordered by marginal likelihood
  expect_true(all(convergence(k3)))
  expect_true(all(diff(marginal_lik(k3)) < 0))
  k2 <- fitModelK(mlist[[3]])
  expect_false(convergence(k2))

  mlist <- list(k3, k2)
  mlist2 <- as(mlist, "MultiBatchList")
  ix <- order(marginal_lik(mlist2), decreasing=TRUE)
  mlist2 <- mlist2[ix]
  expect_identical(modelName(mlist2)[1], "MB3")

  ##k1 <- fitModelK(mlist[[4]])
  burnin(object) <- 50
  iter(object) <- 100
  max_burnin(object) <- 50
  mlist <- mcmc2(object)
  expect_is(mlist, "MultiBatchList")
  expect_identical(modelName(mlist2)[1], "MB3")
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

.test_that("augment data and down sample for MultiBatchList", {

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

test_that("Compute BAF likelihood for each candidate model", {
  skip("uses external data")
  library(SummarizedExperiment)
  mb <- MultiBatchModelExample
  mb <- as(mb, "MultiBatch")
  mapping(mb) <- seq_len(k(mb))
  clist <- CnList(mb)
  data(cnp_se, package="PancCnvsData2")

  hapmap.dir <- system.file("extdata/hapmap", package="PancCnvsData2")
  cnp_57 <- readRDS(file.path(hapmap.dir, "cnp_57.rds"))
  r <- cnp_57[["r"]]
  tidyData <- function(x, value.col){
    probes <- x[[1]]$probeset_id
    x2 <- x %>%
      do.call(cbind, .) %>%
      mutate(probe=probes)
    ix <- c(grep("probeset_id", colnames(x2)),
            grep("plate", colnames(x2)))
    separate <- tidyr::separate
    x3 <- x2[, -ix] %>%
      gather("plate.id", value.col, -probe) %>%
      separate("plate.id", c("plate", "id"), sep="\\.") %>%
      as.tibble
    colnames(x3)[4] <- value.col
    x3
  }
  true.z <- readRDS(file.path(hapmap.dir, "cnp57_truth.rds"))
  truth <- tibble(cn=factor(true.z),
                  id=unique(dat.list[["g"]]$id))
  dat.list <- list(r=tidyData(r, "lrr"),
                   b=tidyData(cnp_57[["b"]], "baf"),
                   g=tidyData(cnp_57[["g"]], "genotype"),
                   truth=truth)
  if(FALSE){
    ## sanity check
    dat <- left_join(dat.list[["b"]], truth, by="id")
    ggplot(dat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  set.seed(132)
  experiment <- hapmapExperiment()
  trace(simulateBatchEffect, browser)
  set.seed(experiment$seed[1])
  hdat <- simulateBatchEffect(dat.list, experiment[1, ])
  snpdat <- hapmapSummarizedExperiment(cnp_57, cnp_se)
  expect_identical(colnames(snpdat), hdat$id)
  if(FALSE){
    ## sanity check
    ggplot(hdat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  colnames(hdat)[1] <- "provisional_batch"
  colnames(hdat)[6] <- "true_batch"
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001)
  snpdat2 <- snpdat[, id(mb)]
  expect_identical(id(mb), colnames(snpdat2))
  mb <- MultiBatchList(data=assays(mb),
                       mp=mcmcParams(mb),
                       iter=100L,
                       burnin=100L,
                       max_burnin=150L)
  mlist <- mcmc2(mb)
  expect_identical(assays(mlist), assays(mb))
  if(FALSE){
    ggMixture(mlist[[1]])
  }
  mb3 <- mlist[[1]]
  expect_identical(modelName(mb3), "MB3")
  clist <- CnList(mb3)
  if(FALSE){
    dat <- tibble(g=genotypes(snpdat)[1, ],
                  b=bafs(snpdat)[1, ],
                  cn=copyNumber(clist[[1]]),
                  true.cn=true.z) %>%
      mutate(cn=factor(cn))
    ggplot(dat, aes(cn, b)) +
      geom_jitter(width=0.1)
    ggplot(dat, aes(true.cn, b)) +
      geom_jitter(width=0.1)
  }
  blik <- modelProb(clist[[1]], snpdat2)
  expect_true(blik > -20)
  stats <- baf_loglik(clist, snpdat2)
  ## model fit looks good, but mapping by baf_loglik does not
  expect_identical(stats$cn.model[1], "0,1,2")
  mb3 <- mlist[[1]]
  mapping(mb3) <- strsplit(stats$cn.model[1], ",")[[1]]
  cn <- factor(copyNumber(mb3))
  pstats <- PancCnvs::performanceStats(assays(mb3)$cn, cn)
})



