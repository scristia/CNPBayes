context("MultiBatch pooled variances")

## constructor
test_that("MultiBatchPooled", {
  model <- MultiBatchModel2()
  expect_is(model, "MultiBatchModel")
  model <- MultiBatchPooled()
  expect_is(model, "MultiBatchPooled")
  expect_equivalent(sigma(model), numeric())

  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -2, -1.95, -0.41, -0.4, -0.395, -0.1,
      0, 0.05), nbatch, k, byrow = FALSE)
  sds <- matrix(0.15, nbatch, k)
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
                            batches=batch(truth))
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
})

.test_that <- function(nm, expr) NULL

.test_that("Marginal likelihood for MultiBatchPooled", {
  data("MultiBatchPooledExample")
  model <- MultiBatchPooledExample
  mp <- McmcParams(iter=500, burnin=1000, nStarts=1, thin=1)
  mcmcParams(model) <- mp
  model2 <- runBurnin(model)
  model2 <- runMcmc(model)
  model2 <- .posteriorSimulation2(model)
  if(FALSE)  ggMultiBatch(model2)

  ## posteriorSimulation returns a sigma matrix
  ##  red_gibbs <- .blockUpdatesBatch(model2, params)
  ##  pstar <- blockUpdates(red_gibbs, root)
  ##ptheta.star <- marginal_theta_batch(model)
  ptheta.star <- theta_multibatch_pvar_red(model2)
  ##model.psigma2 <- reduced_sigma_batch(model.reduced)
  model.psigma2 <- sigma_multibatch_pvar_red(model2)
  expect_true(identical(modes(model.psigma2), modes(model2)))
  ##psigma.star <- p_sigma_reduced_batch(model.psigma2)
  psigma.star <- p_multibatch_pvar_red(model.psigma2)
  ##model.pistar <- reduced_pi_batch(model.reduced)
  model.pistar <- pi_multibatch_pvar_red(model2)
  expect_identical(modes(model.pistar), modes(model2))

  p.pi.star <- p_pmix_reduced_batch(model.pistar)
  ##
  ## Block updates for stage 2 parameters
  ##
  ##model.mustar <- reduced_mu_batch(model.reduced)
  model.mustar <- mu_multibatch_pvar_red(model2)
  expect_identical(modes(model.mustar), modes(model2))
  ## call from multibatch_reduced.cpp
  p.mustar <- p_mu_reduced_batch(model.mustar)
  ##
  ##model.taustar <- reduced_tau_batch(model.reduced)
  model.taustar <- tau_multibatch_pvar_red(model2)
  ##
  expect_identical(modes(model.taustar), modes(model2))
  ## call from multibatch_reduced
  p.taustar <- p_tau_reduced_batch(model.mustar)

  ##model.nu0star <- reduced_nu0_batch(model.reduced)
  model.nu0star <- nu0_multibatch_pvar_red(model2)

  expect_identical(modes(model.nu0star), modes(model2))
  p.nu0star <- prob_nu0_multibatch_pvar_red(model.nu0star)
  ##p.nu0star <- p_nu0_reduced_batch(model.nu0star)

  ##model.s20star <- reduced_s20_batch(model.reduced)
  model.s20star <- s20_multibatch_pvar_red(model2)

  ##p.s20star <- p_s20_reduced_batch(model.s20star)
  p.s20star <- prob_s20_multibatch_pvar_red(model.s20star)

  reduced_gibbs <- cbind(ptheta.star, psigma.star,
                         p.mustar, p.pi.star,
                         p.taustar, p.nu0star,
                         p.s20star)
  colnames(reduced_gibbs) <- c("theta", "sigma", "pi", "mu",
                               "tau", "nu0", "s20")

})

