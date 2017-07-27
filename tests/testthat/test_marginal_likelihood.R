context("Marginal likelihood")

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
  mp <- McmcParams(iter=500, burnin=1000, nStarts=1, thin=5)
  ##    taut2hat = var(theta(model))
  ## qInverseTau2(mn=0.01, sd=0.001)
  mlist <- MarginalModelList(data=galaxies2,
                             mcmc.params=mp,
                             k=1:4,
                             eta.0=200,
                             m2.0=100)
  ## for element of mlist, generate nStarts starting values from a bootstrap
  ## sample and select the best model from the log likelihood evaluated on the
  ## full dataset
  mlist2 <- posteriorSimulation(mlist)
  expect_false(failEffectiveSize(mlist2[[3]]))
  ml <- marginalLikelihood(mlist2)
  expect_equivalent(which.max(ml), 3L)
 
  if(FALSE){
    ## I verified that 1000 burnin, 1000 iter, and thin of 10 converges for k=3
    library(ggplot2)
    p1 <- ggSingleBatchChains(mlist[[2]])
    p1[["comp"]]
    p2 <- ggSingleBatchChains(mlist2[[2]])
    p2[["comp"]]
    p2[["single"]]
    p3 <- ggSingleBatchChains(mlist3[[3]])
    p3[["comp"]]
    library(coda)
    for(i in 1:5){
      ##
      ##
      ##
      mlist2 <- suppressMessages(posteriorSimulation(mlist))
      ##p2 <- ggSingleBatchChains(mlist2[[2]])
      ##p2[["comp"]]
      gd[i, ] <- sapply(mlist2, CNPBayes:::gelmanDiag)
      converged <- gd[i, ] < 1.5
      ##if(converged[3]) stop()
      ml[[i]] <- suppressMessages(marginalLikelihood(mlist2[converged]))

      tmp=.blockUpdates(mlist2[[3]], mcmcParams(mlist2[[3]]))
      tmp <- tmp[tmp[, "theta"] > 0.0001 & tmp[, "sigma"] > 1e-10, ]
      sum(apply(tmp, 2, function(x) log(mean(x))))

      tmp2=.blockUpdates(mlist2[[2]], mcmcParams(mlist2[[2]]))
      tmp2 <- tmp2[tmp2[, "theta"] > 0.0001 & tmp2[, "sigma"] > 1e-10, ]
      sum(apply(tmp2, 2, function(x) log(mean(x))))
    }
  }
})

test_that("batch overfit galaxy", {
  set.seed(1)
  ## load data
  library(MASS)
  data(galaxies)
  ## correct 78th observation
  galaxies[78] <- 26960
  galaxies2 <- (galaxies-median(galaxies))/1000
  galaxies3 <- c(galaxies2, galaxies2 + 10)
  ##mp <- McmcParams(thin=10, iter=1000, burnin=5000, nStarts=1)
  ##mp <- McmcParams(thin=10, iter=1000, burnin=50, nStarts=100)
  mp <- McmcParams(burnin=200, nStarts=20, iter=100)
  model.list <- BatchModelList(data=galaxies3,
                               batch=rep(1:2, each=length(galaxies)),
                               k=1:4,
                               mcmc.params=mp,
                               eta.0=0.08,
                               m2.0=50)
  ## default prior on tau is far too informative for the galaxy data
  ##hypp <- HyperparametersBatch(eta.0=0.08, m2.0=50, k=3)
  mlist <- posteriorSimulation(model.list)
  pstar <- marginal_theta_batch(mlist[[3]])
  expect_false(failSmallPstar(pstar))
  pstar4 <- marginal_theta_batch(mlist[[4]])
  ##pstar <- .blockUpdatesBatch(mlist[[3]], mlParams())
  ##expect_true(failSmallPstar(pstar4))
  ##expect_warning(.blockUpdatesBatch(mlist[[4]], mlParams()))
  ##sum(round(log(apply(pstar, 2, mean)), 3))
  ml <- marginalLikelihood(mlist, mlParams(ignore.effective.size=TRUE,
                                           warnings=FALSE))
  expect_equivalent(which.max(ml), 3L)
  ## For the k=4 model, there is some label switching. A correction factor is
  ## not needed in the calculation of the marginal likelihood.
  expect_false(failEffectiveSize(mlist[[3]]))
  if(FALSE){
    ## verify stage two log lik is reasonable for the precision
    ## (code below from overview vignette)
    mod <- mlist.mb[[2]]
    mod <- useModes(mod)
    prec <- 1/sigma2(mod)
    shape <- 0.5*nu.0(mod)
    rate <- 0.5*nu.0(mod) * sigma2.0(mod)
    dgamma(prec, shape=shape, rate=rate)
    ## doesn't make sense that modal value of prec would have near zero prob.
    ## with modal values of nu.0 and sigma2.0...
    hist(rgamma(1000, shape=shape, rate=rate), breaks=100)
    ##likprec[0] = sum(log(dgamma(1.0/sigma2(_, k),
    ##                            0.5*nu0[0], 1.0/(0.5*nu0[0]*s20[0])))) ;
    plist <- ggMultiBatchChains(mlist[[3]])
    plist[["batch"]]
    plist4 <- ggMultiBatchChains(mlist[[4]])
    plist4[["batch"]]
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
  model <- MarginalModel(data=galaxies2,
                         hypp=hypp,
                         mcmc.params=mp)


  nstarts <- 100
  mlist <- replicate(100, MarginalModel(y(model),
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
