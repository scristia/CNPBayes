test_relabel_marginal <- function(){
  set.seed(1)
  truth <- simulateData(N=2500, p=rep(1/3, 3),
                        theta=c(-1, 0, 1),
                        sds=rep(0.1, 3))
  mp <- McmcParams(iter=500, burnin=500)
  model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
  model <- posteriorSimulation(model)

  ## relabel z:  3->1 and 1->3
  zindex <- c(3, 2, 1)
  model2 <- relabel(model, zindex)
  checkIdentical(as.integer(tablez(model)), as.integer(tablez(model2)[zindex]))
  checkIdentical(CNPBayes:::dataMean(model), CNPBayes:::dataMean(model2)[zindex])
  ##
  ## All other parameters should be the same
  ##
  checkIdentical(theta(model), theta(model2))
  checkIdentical(mu(model), mu(model2))
  checkIdentical(p(model), p(model2))
  checkIdentical(sigma(model), sigma(model2))
  checkIdentical(nu.0(model), nu.0(model2))
  if(FALSE){
    ##
    ## should observe chain swap
    ##
    burnin(model2) <- 0
    model2 <- posteriorSimulation(model2)
    head(thetac(model2))
    plot.ts(thetac(model2), plot.type="single", col=1:3)
    set.seed(1)
    truth <- simulateData(N=2500, p=rep(1/3, 3),
                          theta=c(-1, 0, 1),
                          sds=rep(0.5, 3))
    mp <- McmcParams(iter=500, burnin=500)
    model <- MarginalModel(data=y(truth), k=3, mcmc.params=mp)
    model <- posteriorSimulation(model)
    zindex <- c(3, 2, 1)
    model2 <- relabel(model, zindex)
    burnin(model2) <- 0
    model2 <- posteriorSimulation(model2)
    par(mfrow=c(1,2))
    plot.ts(thetac(model), col=1:3, plot.type="single")
    plot.ts(thetac(model2), col=1:3, plot.type="single")
  }
}
