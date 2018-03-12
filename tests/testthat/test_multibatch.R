context("MultiBatchModel")

.test_that <- function(name, expr)  NULL

test_that("initial values", {
  extdata <- system.file("extdata", package="CNPBayes")
  dat <- readRDS(file.path(extdata, "mckean_data.rds"))
  data <- dat$y
  batch <- dat$batch
  mcmc.params <- McmcParams()
  k <- 1:4
  model.list <- vector("list", length(k))
  set.seed(2)
  j <- 3
  hypp <- HyperparametersMultiBatch(k=4)
  model <- MB(dat=dat$y,
              batches=as.integer(factor(dat$batch)),
              hp=hypp)
  expect_identical(ncol(sigma(model)), 4L)
  expect_true(!is.na(log_lik(model)))
})

.test_that("MultiBatchModel", {
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
      0, 0.05), nbatch, k, byrow = FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[, 1] <- 0.3
  N <- 1000
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  batches <- batch(truth)
  dat <- y(truth)
  hp <- HyperparametersMultiBatch(k=3,
                             mu=-0.75,
                             tau2.0=0.4,
                             eta.0=32,
                             m2.0=0.5)
  ##trace(MultiBatchModel, browser)
  expect_true(validObject(MB()))
  mp <- McmcParams(iter = 1000,
                   burnin = 1000,
                   nStarts = 4,
                   thin=10)
  model <- MB(hp=hp, mp=mp, dat=y(truth),
                           batches=batch(truth))
  expect_true(validObject(model))

  library(purrr)
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts = 1, thin=1)
  m <- MB(dat=y(truth),
                        mp=mp, hp=hp,
                        batches=batch(truth))
  m2 <- posteriorSimulation(m)
  gelman_rubin(mcmcList(list(m2)), hp)

  mod.list <- replicate(4, MB(dat=y(truth),
                                           mp=mp, hp=hp,
                                           batches=batch(truth)))
  mod.list2 <- map(mod.list, posteriorSimulation)
  model <- combine_batch(mod.list2, batches=batch(mod.list2[[1]]))
  mp <- McmcParams(iter = 1000, burnin = 1000, nStarts = 4, thin=5)
  set.seed(4894)
  model <- gibbs_batch(dat=y(truth), mp=mp, hp=hp,
                       batches=batch(truth))
  trace(r_theta_multibatch, browser)
  r_theta_multibatch(model)
  compute_marginal_lik(model)
  .ml_batchmodel(model)

  marginal_lik(model)

  k(hp) <- 4
  m4 <- gibbs_batch(dat=y(truth), mp=mp, hp=hp,
                    batches=batch(truth))
  expect_true(is.na(marginal_lik(m4)))

  k(hp) <- 1
  mod <- MB(dat=y(truth), batches=batch(truth), hp=hp, mp=mp)
  k(hp) <- 2
  mod <- MB(dat=y(truth), batches=batch(truth), hp=hp, mp=mp)
  mod.list <- gibbs_batch_K(dat=y(truth),
                            mp=mp, hp=hp,
                            k_range=c(1, 4),
                            batches=batch(truth))

  hp.sb <- Hyperparameters(tau2.0=0.4,
                           mu.0=-0.75,
                           eta.0=32,
                           m2.0=0.5)
  ## ignoring the batches
  mod.list2 <- gibbs_K(dat=y(truth), mp=mp, hp=hp.sb)
  hp.list <- list(single_batch=hp.sb,
                  multi_batch=hp)
  models <- gibbs_all(hp.list=hp.list, dat=y(truth),
                      batches=batch(truth),
                      mp=mp,
                      top=3)
  expect_is(models[[1]], "MultiBatchModel")
  expect_identical(k(models[[1]]), 3L)
  p1 <- ggMultiBatch(truth)
  p2 <- ggMultiBatch(models[[1]])
  library(gridExtra)
  grid.arrange(p1, p2)
  ## ML_MB2 > ML_SB3
  ggMultiBatch(models[[1]])
  ggMultiBatch(models[[2]])
  ggSingleBatch(models[[3]])
  purrr::map_dbl(models, marginal_lik)
})


test_that("easy", {
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
      nbatch, k, byrow = FALSE)
  sds <- matrix(0.1, nbatch, k)
  N <- 1500
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             theta = means, sds = sds,
                             p = c(1/5, 1/3, 1 - 1/3 - 1/5))
  ##yy <- y(truth)
  ##expect_identical(yy[order(batch(truth))], yy)
  mcmcp <- McmcParams(iter = 50, burnin = 0)
  set.seed(123)
  model <- MB(dat=y(truth), batches = batch(truth),
              hp=hpList(k = 3)[["MB"]],
              mp = mcmcp)
  u1 <- u(model)
  model <- posteriorSimulation(model)
  u2 <- u(model)
  expect_true(!identical(u1, u2))

  model <- startAtTrueValues(model, truth)
  expect_identical(batch(truth), batch(model))
  expect_identical(y(truth), y(model))
  expect_identical(theta(truth), theta(model))
  expect_identical(sigma(truth), sigma(model))
  expect_identical(p(truth), p(model))
  expect_identical(z(truth), z(model))
  iter(model) <- 50L
  if (FALSE) {
      model <- .Call("mcmc_batch", model, mcmcParams(model))
      set.seed(123)
      .Call("update_sigma20_batch", model)
      set.seed(123)
      .updateSigma2.0Batch(model2)
      eta.0 <- 1800/100
      m2.0 <- 1/(60/100)
      a <- 0.5 * eta.0
      b <- 0.5 * eta.0 * m2.0
      x <- rgamma(1000, a, rate = 1/b)
      ix <- 1/x
      hist(sqrt(1/x), breaks = 100)
      s <- mean(sqrt(1/x))
      theta <- rnorm(1000, -1, s)
      hist(theta, breaks = 100, xlim = c(-2, 1.5))
  }
  set.seed(1)
  mcmcp <- McmcParams(iter = 300, burnin = 300, nStarts = 5)
  model <- MB(dat=y(truth),
                            batches = batch(truth),
                            hp=hpList(k=3)[["MB"]],
                            mp=mcmcp)
  model <- posteriorSimulation(model)
  expect_equal(theta(model), theta(truth),
               scale=0.01, tolerance=1)
  if (FALSE) {
    MultiBatchModelExample <- model
    save(MultiBatchModelExample, file = "data/MultiBatchModelExample.RData")
  }
  ##
  ## test upSample2 without upsampleing
  ##
  ##  - z-probabilities obtained from the theoretical distribution should be close to the actual posterior probabilties
  ##
  tab <- tibble(medians=y(model), batch_index=batch(model))
  model.up <- upSample2(orig.data=tab, model, up_sample=FALSE)
  pz.up <- probz(model.up)
  pz <- probz(model)
  expect_equal(pz, pz.up, tolerance=0.01, scale=1)
})




hardTruth <- function(p1, s, N=1000){
  set.seed(1234)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.9, -2, -1.85,
                    -0.45, -0.4, -0.35,
                    -0.1, 0, -0.05), nbatch, k, byrow=FALSE)
  sds <- matrix(s, nbatch, k)
  ##p1 <- prop_comp1
  p2 <- 20*p1
  p3 <- 1-p1-p2
  ##N <- 2000
  truth <- simulateBatchData(N,
                             theta=means,
                             sds=sds,
                             batch=rep(letters[1:3], length.out=N),
                             p=c(p1, p2, p3))
  truth
}

test_that("stay_near_truth", {
  set.seed(123)
  library(SummarizedExperiment)
  ## embed function in test for now
  # end function
  truth <- hardTruth(p1=0.02, s = 0.1)
  se <- as(truth, "SummarizedExperiment")
  ##
  ## nStarts=0 forces chain to start at true values
  ##
  ## - these unit tests verify that the model stays in a region of high
  ## - posterior prob.
  mcmcp <- McmcParams(iter = 100, burnin = 0, nStarts=0)
  modelk <- MB(dat = y(truth), batches = batch(truth),
                             hp=hpList(k = 3)[["MB"]],
                             mp = mcmcp)
  modelk <- startAtTrueValues(modelk, truth)
  mmodel <- posteriorSimulation(modelk)
  pmns <- thetaMean(mmodel)
  ##ps <- CNPBayes:::sigmaMean(mmodel)
  s <- sigma(mmodel)
  pmix <- pMean(mmodel)
  expect_equal(theta(truth), pmns, tolerance=0.04)
  expect_equal(p(truth), pmix, tolerance=0.04)
  ## TODO: This example could be extended to focus on what to do when a batch
  ## does not have any homozygous deletions. Batch 2 has only 1 homozygous
  ## deletion
})

test_that("kbatch", {
  set.seed(123)
  k <- 3
  means <- matrix(c(rnorm(5, -1, 0.1), rnorm(5, 0, 0.1), rnorm(5,
      1, 0.1)), 5, k, byrow = FALSE)
  sds <- matrix(0.1, 5, k)
  N <- 3000
  ## the last 2 batches are much smaller
  probs <- c(1/3, 1/3, 3/10, 0.02, 0.013)
  probs <- probs/sum(probs)
  batch <- sample(1:5, size = N, prob = probs, replace = TRUE)
  p <- c(1/5, 1/3)
  p <- c(p, 1 - sum(p))
  ##trace(simulateBatchData, browser)
  truth <- simulateBatchData(N = N, batch = batch, theta = means,
      sds = sds, p = p)
  mp <- McmcParams(iter = 100, burnin = 50, nStarts = 10)
  kmod <- MB(dat=y(truth),
             batches=batch(truth),
             hp=hpList(k = 3)[["MB"]],
             mp = mp)
  kmod <- posteriorSimulation(kmod)
  expected <- max.col(probz(kmod))
  cn <- map_z(kmod)
  expect_identical(cn, expected)

  set.seed(1000)
  index <- sort(unique(c(sample(seq_len(N), 500), which(batch(kmod) %in% 4:5))))
  mp <- McmcParams(iter = 100, burnin = 100, nStarts = 10)
  kmod2 <- MB(dat=y(kmod)[index],
              batches=batch(kmod)[index],
              hp=hpList(k=3)[["MB"]],
              mp=mp)
  mp <- McmcParams(iter = 100, burnin = 0, nStarts = 10)
  kmod3 <- posteriorSimulation(kmod2)
  cn2 <- map_z(kmod3)
  pz <- probz(kmod3)
  pz <- mapCnProbability(kmod3)
  tab.z <- as.integer(table(z(kmod3)))
  tab.z2 <- colSums(round(pz, 1))
  expect_equal(tab.z, tab.z2, tolerance=1)
  if (FALSE) {
      fit <- list(posteriorSimulation(kmod, k = 1), posteriorSimulation(kmod,
          k = 2), posteriorSimulation(kmod, k = 3), posteriorSimulation(kmod,
          k = 4))
      fit <- marginalLikelihood(fit)
      prz <- probz(fit$models[[4]])
      cn <- map_z(fit$models[[4]])
      plot(r, cn, pch = 20, cex = 0.3)
      prz <- cnProbability(prz, 4)
      plot(jitter(prz, amount = 0.05), jitter(cn, amount = 0.05),
           pch = 20, cex = 0.3)
      table(cn)
      pz <- cnProbability(probz(fit$models[[4]]), 4)
      r <- y(fit$models[[4]])
      plot(r, pz, pch = ".")
      expect_true(k(orderModels(fit))[1] == 3)
  }
})

test_that("test_unequal_batch_data", {
    expect_error(MB(dat = 1:10, batches = 1:9))
})
