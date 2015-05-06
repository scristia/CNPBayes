test_loglik <- function(){
  set.seed(2000)
  library(oligoClasses)
  truth <- simulateData(N=5e3,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  ll.truth <- logLik(truth)
  checkEquals(ll.truth, computeLoglik(truth))

  yy <- y(truth)
  th <- theta(truth)
  sd <- sigma(truth)
  p_ <- p(truth)
  log_lik <- sum(log(p_[1]*dnorm(yy, th[1], sd[1]) + p_[2]*dnorm(yy, th[2], sd[2]) +
                         p_[3]*dnorm(yy, th[3], sd[3])))
  checkEquals(log_lik, ll.truth)

  ##
  ##  Simulated batch variable has nothing to do with the data.
  ##
  model <- BatchModel(data=y(truth),
                      batch=rep(letters[1:3], length.out=length(y(truth))),
                      k=3)
  iter(model, force=TRUE) <- 250
  burnin(model) <- 100
  nStarts(model) <- 5
  model <- posteriorSimulation(model)
  ll1 <- logLik(model)
  ll2 <- computeLoglik(model)
  checkEquals(ll1, ll2, tolerance=0.005)

  yy <- y(model)
  th <- theta(model)
  sd <- sigma(model)
  p_ <- table(batch(model), z(model))
  p_ <- p_/rowSums(p_)
  ##p_ <- p(model)
  ub <- unique(batch(model))
  log_lik <- rep(0, length(y(model)))
  for(b in ub){
    this_batch <- batch(model) == ub[b]
    log_lik <- log_lik + log(p_[b, 1]*dnorm(yy, th[b, 1], sd[b, 1]) +
                                 p_[b, 2]*dnorm(yy, th[b, 2], sd[b, 2]) +
                                     p_[b, 3]*dnorm(yy, th[b, 3], sd[b, 3])) * this_batch
  }
  ll3 <- sum(log_lik)
  checkEquals(ll2, ll3)
}

test_senseless_batch <- function(){
  set.seed(2000)
  library(oligoClasses)
  truth <- simulateData(N=2500,
                        theta=c(-2, -0.4, 0),
                        sds=c(0.3, 0.15, 0.15),
                        p=c(0.005, 1/10, 1-0.005-1/10))
  set.seed(123)
  m1 <- computeMarginalLik(y(truth))
  models <- orderModels(m1)
  checkTrue(k(models[[1]]) == 3)
  set.seed(123)

  ##
  ## TODO: fix batch model
  ##
  mp <- McmcParams(iter=10, burnin=10, nStarts=1)
  b <- marginal(object=se, batch=batch(truth), mcmc.params=mp, maxperm=3)
  checkTrue(is(b, "MarginalModelList"))
  mp <- McmcParams(iter=300, burnin=300, nStarts=10)
  senseless_batch <- sample(1:3, length(y(truth)), replace=TRUE)
  b <- marginal(object=se, batch=senseless_batch, mcmc.params=mp, maxperm=3)
  bf <- bayesFactor(rbind(summary(m1), summary(b)))
  checkTrue(substr(names(bf), 1, 2) == "M3")
}
