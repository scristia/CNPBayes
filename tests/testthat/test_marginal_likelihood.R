context("Marginal likelihood")

test_that("overfit model", {
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
  mlist <- posteriorSimulation(model, k=2:4)
  mcmcParams(mlist) <- McmcParams(nStarts=0, burnin=2000, iter=1000, thin=0)
  mlist2 <- posteriorSimulation(mlist)
  ml <- marginalLikelihood(mlist2)
  ## force calculation of marginal likelihood for overfit model
  ml2 <- marginalLikelihood(mlist2, mlParams(ignore.effective.size=TRUE,
                                             warnings=FALSE))
  expect_equivalent(which.max(ml2), 2)
  ##gd <- matrix(NA, 5, length(mlist))
  if(FALSE){
    ## I verified that 1000 burnin, 1000 iter, and thin of 10 converges for k=3
    library(ggplot2)
    p1 <- ggSingleBatchChains(mlist[[2]])
    p1[["comp"]]
    p2 <- ggSingleBatchChains(mlist2[[2]])
    p2[["comp"]]
    p2[["single"]]
    p3 <- ggSingleBatchChains(mlist2[[3]])
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
  mp <- McmcParams(burnin=0, nStarts=1000)
  ## default prior on tau is far too informative for the galaxy data
  hypp <- HyperparametersBatch(eta.0=0.08, m2.0=50, k=3)
  model <- BatchModel(data=galaxies3,
                      batch=rep(1:2, each=length(galaxies)),
                      mcmc.params=mp, k=3,
                      hypp=hypp)
  mlist <- posteriorSimulation(model, k=2:4)
  mcmcParams(mlist) <- McmcParams(thin=5, nStarts=0, iter=1000, burnin=0)
  mlist2 <- posteriorSimulation(mlist)
  mlist3 <- posteriorSimulation(mlist2)
  ##
  ## visual inspection of k=3 for mlist3 shows that their is no label swapping
  ## -- expect p(theta* | ...) >> 0
  ##
  ##plist3 <- ggMultiBatchChains(mlist3[[2]])
  ##plist3[["batch"]]
  pstar <- marginal_theta_batch(mlist3[[2]])
  expect_false(failSmallPstar(pstar))
  pstar4 <- marginal_theta_batch(mlist3[[3]])


  pstar <- .blockUpdatesBatch(mlist3[[2]], mlParams())
  expect_warning(.blockUpdatesBatch(mlist3[[3]], mlParams()))
  sum(round(log(apply(pstar, 2, mean)), 3))
  ml <- marginalLikelihood(mlist3, mlParams(ignore.effective.size=TRUE,
                                            warnings=FALSE))
  expect_equivalent(which.max(ml), 2L)
  ## For the k=4 model, there is some label switching. A correction factor is
  ## not needed in the calculation of the marginal likelihood.
  expect_true(failEffectiveSize(mlist3[[3]]))
  expect_false(failEffectiveSize(mlist3[[2]]))

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
    plist <- ggMultiBatchChains(mod)
    plist[["single"]]
  }
})
