context("MultiBatch pooled variances")

## constructor
test_that("MultiBatchPooled", {
  model <- MultiBatchModel2()
  expect_is(model, "MultiBatchModel")
  model <- MultiBatchPooled()
  expect_is(model, "MultiBatchPooled")
  expect_equivalent(sigma(model), numeric())
  mp <- McmcParams(iter=50, burnin=5)
  model <- MultiBatchPooledExample
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  model2 <- as(model, "MultiBatchCopyNumberPooled")
  expect_true(validObject(model2))
  model3 <- useModes(model2)

  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
                    0, 0.05),
                  nbatch, k, byrow = FALSE)
  sds <- matrix(rep(c(0.1, 0.15, 0.2), each=3),
                nrow=nbatch, ncol=k, byrow=TRUE)
  N <- 300
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  if(FALSE)  ggMultiBatch(truth)

  hp <- HyperparametersMultiBatch(k=3,
                                  mu=-0.75,
                                  tau2.0=0.4,
                                  eta.0=32,
                                  m2.0=0.5)
  model <- MultiBatchPooled(dat=y(truth),
                            hp=hp,
                            batches=batch(truth),
                            mp=mp)

  u1 <- u(model)
  model <- posteriorSimulation(model)
  u2 <- u(model)
  expect_true(!identical(u1, u2))

  expect_true(validObject(model))
  if(FALSE){
    MultiBatchPooledExample <- model
    save(MultiBatchPooledExample, file="data/MultiBatchPooledExample.rda")
  }
  expect_is(sigma(model), "numeric")
  expect_identical(length(sigma(model)), 3L)
})

test_that("MultiBatchPooled MCMC", {
  data("MultiBatchPooledExample")
  model <- MultiBatchPooledExample
  pr <- multinomialPr_multibatch_pvar(model)
  zz <- z_multibatch_pvar(model)
  z(model) <- zz
  ## unchanged
  mns <- compute_means_batch(model)
  dataMean(model) <- mns
  prec <- compute_prec_batch(model)
  dataPrec(model) <- prec
  thetas <- theta_multibatch_pvar(model)
  theta(model) <- thetas
  sigma2s <- sigma2_multibatch_pvar(model)
  sigma2(model) <- sigma2s
  ## same as MultiBatch
  mus <- update_mu_batch(model)
  mu(model) <- mus
  ## unchanged
  tau2s <- update_tau2_batch(model)
  tau2(model) <- tau2s

  s20 <- sigma20_multibatch_pvar(model)
  sigma2.0(model) <- s20
  nu0 <- nu0_multibatch_pvar(model)
  nu.0(model) <- nu0
  ## unchanged
  ps <- update_p_batch(model)

  ll <- loglik_multibatch_pvar(model)
  ll <- stagetwo_multibatch_pvar(model)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=2)
  model2 <- burnin_multibatch_pvar(model, mp)
  model2 <- mcmc_multibatch_pvar(model, mp)
  expect_true(validObject(model2))
})

.test_that <- function(nm, expr) NULL


test_that("Marginal likelihood for MultiBatchPooled", {
  library(magrittr)
  library(tibble)
  data("MultiBatchPooledExample")
  model <- MultiBatchPooledExample
  set.seed(123)
  mp <- McmcParams(iter=50, burnin=10, nStarts=1, thin=1)
  mcmcParams(model) <- mp
  model2 <- runBurnin(model)
  model2 <- runMcmc(model2)
  model2 <- .posteriorSimulation2(model2)
  if(FALSE)  ggMultiBatch(model2)

  ptheta.star <- theta_multibatch_pvar_red(model2)
  model.psigma2 <- sigma_multibatch_pvar_red(model2)
  expect_true(identical(modes(model.psigma2), modes(model2)))
  psigma.star <- psigma_multibatch_pvar_red(model.psigma2)
  model.pistar <- pi_multibatch_pvar_red(model2)
  expect_identical(modes(model.pistar), modes(model2))
  p.pi.star <- p_pmix_reduced_batch(model.pistar)
  ##
  ## Block updates for stage 2 parameters
  ##
  model.mustar <- mu_multibatch_pvar_red(model2)
  expect_identical(modes(model.mustar), modes(model2))
  p.mustar <- p_mu_reduced_batch(model.mustar)
  model.taustar <- tau_multibatch_pvar_red(model2)
  expect_identical(modes(model.taustar), modes(model2))
  p.taustar <- p_tau_reduced_batch(model.mustar)
  model.nu0star <- nu0_multibatch_pvar_red(model2)
  expect_identical(modes(model.nu0star), modes(model2))
  p.nu0star <- pnu0_multibatch_pvar_red(model.nu0star)
  model.s20star <- s20_multibatch_pvar_red(model2)
  p.s20star <- ps20_multibatch_pvar_red(model.s20star)
  reduced_gibbs <- cbind(ptheta.star, psigma.star,
                         p.mustar, p.pi.star,
                         p.taustar, p.nu0star,
                         p.s20star)
  colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                               "tau", "nu0", "s20")

  set.seed(15)
  model3 <- runBurnin(model)
  model3 <- runMcmc(model3)
  model3 <- .posteriorSimulation2(model3)
  params <- mlParams(ignore.small.pstar=TRUE)
  pmat <- .blockUpdatesMultiBatchPooled(model3, params) %>%
    as.tibble
  ml <- .ml_multibatch_pooled(model3, params)
  ##expect_equal(as.numeric(ml), -57, tolerance=0.05)
  expect_equal(as.numeric(ml), -108, tolerance=0.05)
  marginal_lik(model3) <- ml
})

test_that("MultiBatchPooled model selection", {
  set.seed(123)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
                    0, 0.05),
                  nbatch, k, byrow = FALSE)
  sds <- matrix(rep(c(0.1, 0.15, 0.2), each=3),
                nrow=nbatch, ncol=k, byrow=TRUE)
  N <- 500
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  mp <- McmcParams(iter=1000, burnin=1000, nStarts=4, thin=1)
  hp <- HyperparametersMultiBatch(k=4,
                                  mu=-0.75,
                                  tau2.0=0.4,
                                  eta.0=32,
                                  m2.0=0.5)
  batches <- batch(truth)
  set.seed(941)
  ##options(warn=2, error=utils::recover)
  model <- MultiBatchPooled(dat=y(truth), mp=mp, hp=hp,
                            batches=batch(truth))
  ## fit model with k=4
  ##
  ## running too few iterations for this to be very useful
  expect_warning(
    model <- gibbs_multibatch_pooled(hp,
                                     mp=McmcParams(iter=100, burnin=100, nStart=4),
                                     y(truth),
                                     batches=batch(truth))
  )
  if(FALSE){
    figs <- ggChains(model)
    figs$theta
    figs$sigma
    ## why is there so much variation in the mus?
    figs$comp
    figs$sigma
    ggMultiBatch(model)
  }
  ## fit models k=1 -> k=4
  if(FALSE){
    mlist <- gibbsMultiBatchPooled(hp=hp,
                                   mp=mp,
                                   k_range=c(3, 3),
                                   dat=y(truth),
                                   batches=batch(truth))
    model <- mlist[[1]]
    ggMixture(model)
    tab <- posteriorPredictive(model)
    ggPredictive(model, tab)
    expect_identical(k(model), 3L)
  }
  if(FALSE){
    mlist <- gibbs("MBP",
                   hp.list=list(MBP=hp),
                   mp=mp,
                   k_range=c(3, 3),
                   dat=y(truth),
                   batches=batch(truth))
    model <- mlist[[1]]
    ggPredictive(model)
    ggMixture(model)
    figs <- ggChains(model)
    figs$theta
    figs$sigma
    figs$comp
    figs$sigma
  }
})
