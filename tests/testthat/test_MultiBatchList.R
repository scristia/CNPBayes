contect("MultiBatch classes")

test_that("revised_constructors", {
  mb <- MultiBatchModel2()
  expect_true(validObject(mb))

  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mb2 <- as(sb, "MultiBatch")

  ##
  ## Initialize from data
  ##
  tab <- assays(mb2)
  tmp <- MultiBatch(data=tab)


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
  ## TODO chains as a list

##  S <- 10; K <- 3; B <- 1
##  mcmc.list(matrix(NA, S, K*B),
##            matrix(NA, S, K*B)),
##  pi=matrix(NA, S, K),
##  mu=numeric(S),
##  tau2=numeric(S),
##  nu.0=numeric(S),
##  sigma2.0=numeric(S),
##  logprior=numeric(S),
##  loglik=numeric(S),
##  zfreq=matrix(as.integer(NA), S, K))

  ##
  ## MultiBatchList
  ## @assays
  ## @summaries.list  ## each element corresponds to k
  ## @chains.list
  ## @parameters.list
  ## @flags.list
  ## converged(MBL) ##vector
  ## marginal_lik(MBL) ##vector
  ## sort(MBL) ## arranges list elements by marginal like
  ## MBL[ converged(MBL) ] ## only converged models
  ##
  ## MB <- MBL[[1]] ##  MultiBatch model
  ## modelNames(MBL)
  ## MB <- MBL[["MB"]]
  
})
