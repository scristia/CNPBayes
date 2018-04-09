context("Marginal likelihood")

.test_that <- function(nm, expr) NULL

test_that("overfit model", {
  ## dataset is small  -- priors are more influential
  ## -- default priors set for CNP data are less effective for the galaxy data
  ##    namely, the prior on tau2
  ##    taut2hat = var(theta(model))
  set.seed(1)
  library(MASS)
  data(galaxies)
  # correct 78th observation
  galaxies[78] <- 26960
  galaxies2 <- (galaxies-median(galaxies))/1000
  ##
  ## stronger prior on variance is needed to reduce label switching
  ##
  mp <- McmcParams(iter=1000, burnin=1000, nStarts=5, thin=2)
  hp <- Hyperparameters(k=2,
                        mu=-0.75,
                        tau2.0=0.4,
                        eta.0=200, ## default 32
                        m2.0=100)  ## default 0.5
  if(FALSE){
    model <- gibbs_K(hp=hp, mp=mp, dat=galaxies2, k_range=c(1, 3),
                     max_burnin=1000)
    ##    taut2hat = var(theta(model))
    ## qInverseTau2(mn=0.01, sd=0.001)
    ggChains(model[[1]])
    ggSingleBatch(model[[1]])
    expect_identical(names(model)[1], "SB3")
  }
})

.test_that <- function(name, expr) NULL
.test_that("number starts", {
  set.seed(1)
  # load data
  library(MASS)
  data(galaxies)
  # correct 78th observation
  galaxies[78] <- 26960
  mp <- McmcParams(iter=0, burnin=0, nStarts=20)
  ##
  ## Need to adjust hyper-parameters for the galaxy data
  ##
  galaxies2 <- (galaxies-median(galaxies))/sd(galaxies)
  hypp <- Hyperparameters(type="marginal", eta.0=200, m2.0=50,
                          a=2, b=1, tau2.0=1000)
  model <- SingleBatchModel(data=galaxies2,
                         hypp=hypp,
                         mcmc.params=mp)


  nstarts <- 100
  mlist <- replicate(100, SingleBatchModel(y(model),
                                        mcmc.params=mcmcParams(model),
                                        hypp=hyperParams(model),
                                        k=k(model)))

  thetas <- lapply(mlist, theta)
  thetas <- do.call(cbind, thetas)



  mlist <- posteriorSimulation(model, k=2:4)
  mcmcParams(mlist) <- McmcParams(nStarts=0, burnin=2000, iter=1000, thin=0)
  mlist2 <- posteriorSimulation(mlist)
  ml <- marginalLikelihood(mlist2)
  ## force calculation of marginal likelihood for overfit model
  ml2 <- marginalLikelihood(mlist2, mlParams(ignore.effective.size=TRUE,
                                             warnings=FALSE))
  expect_equivalent(which.max(ml2), 2)
})
