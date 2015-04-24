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
  iter(model) <- 250
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

  se <- as(truth, "SummarizedExperiment")
  mp <- mcmcParams(model)
  iter(mp) <- 250
  nStarts(mp) <- 10
  set.seed(123)
  m1 <- marginal(se, mcmc.params=mp, maxperm=3)
  if(FALSE){
    plot.ts(thetac(m1[[4]]), plot.type="single", col=1:4)
  }
  iter(m1) <- c(0, 250, 300, 300)
  burnin(m1) <- c(0, 250, 350, 350)
  nStarts(m1) <- 1
  m1 <- marginal(m1)
  my <- summary(m1)
  checkTrue(best(m1) == "3")
  set.seed(123)
  mp <- McmcParams(iter=300, burnin=300, nStarts=10)
  m2 <- marginal(se, batch=batch(model), mcmc.params=mp, maxperm=3)
  if(FALSE){
    ##pc <- prcomp(yy, center=TRUE, scale.=TRUE)$x[, 1]
    ##if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
    plot(m2[[3]])
    par(mfrow=c(1, 3), las=1)
    tracePlot(m2[[3]], "theta", col=1:3)
    tracePlot(kmod, "theta", col=1:3, ylim=c(-3,2))
    tmp <- thetac(m2[[3]])
  }
  by <- summary(m2)
  iter(m2) <- c(0, 200, 500, 500)
  burnin(m2) <- c(0, 200, 500, 500)
  nStarts(m2) <- 1
  load_all()
  m2 <- marginal(m2)
  by <- summary(m2)
  ## Check that 4 component model is unstable
  checkTrue(by["B4", "range"] > 10)
  ## Check marginal for 3 component model is consistent
  checkTrue(by["B3", "range"] < 7)
  my <- rbind(my, by)
  ##trace(bayesFactor, browser)
  bf <- bayesFactor(my)
  ## Select the marginal model with 3 components over the batch model
  ## with 3 components
  checkIdentical(names(bf), "M3-B3")
}
