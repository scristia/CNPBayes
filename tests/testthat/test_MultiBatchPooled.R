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
  ## I think this could be replaced by the method for MultiBatchModel
  mns <- means_multibatch_pvar(model)
  dataMean(model) <- mns
  ##
  ## SingleBatchPooled used compute_prec defined for SingleBatchModel
  ## -- the precision is a k-length vector
  ## -> left vars_multibatch_pvar the same as in MultiBatch
  ## -> left prec_multibatch_pvar the same as in MultiBatch
  vars <- vars_multibatch_pvar(model)
  prec <- prec_multibatch_pvar(model)
  dataPrec(model) <- prec
  thetas <- theta_multibatch_pvar(model)
  theta(model) <- thetas
  sigma2s <- sigma2_multibatch_pvar(model)
  sigma2(model) <- sigma2s
  mus <- mu_multibatch_pvar(model)
  tau2s <- tau2_multibatch_pvar(model)
  s20 <- sigma20_multibatch_pvar(model)
  nu.0 <- nu0_multibatch_pvar(model)
  ps <- p_multibatch_pvar(model)
})
