context("Marginal likelihood calibration")

.test_that <- function(nm, expr) NULL

.test_that("galaxy model", {
    # set seed
    set.seed(1)
    # load data
    library(MASS)
    data(galaxies)
    # correct 78th observation
    galaxies[78] <- 26960
    # fit model
    mp <- McmcParams(thin=10, iter=1000, burnin=10000, nStarts=5)
    hp <- Hyperparameters(k=3, m2.0=100)
    ## specified burnin is more than the max_burnin
    expect_error(gibbs_K("SB",
                         mp=mp, hp=hpList(k=3)[["SB"]],
                         k_range=c(3, 3), dat=galaxies/1000,
                         max_burnin=10000))
    model <- gibbs("SB",
                   mp=mp,
                   hpList(k=3),
                   k_range=c(3, 3),
                   dat=galaxies/1000,
                   max_burnin=burnin(mp))[[1]]
    ggMixture(model)
    ggChains(model)[[1]]
    ml <- marginal_lik(model)
    # calculate marginal likelihood and compare to "truth"
    published.mlik <- -226.791
    m.y <- unname(ml)
    expect_equal(m.y, published.mlik, tolerance=3, scale=1)
})


.test_that("scratch", {
  nchains <- nStarts(mp)
  if(nchains==1) stop("Must initialize at least 2 chains with nStarts ")
  nStarts(mp) <- 1L ## because posteriorsimulation uses nStarts in a different way
  if(iter(mp) < 500){
    stop("Require at least 500 Monte Carlo simulations")
  }
  if(burnin(mp) > max_burnin) stop("Specified burnin is greater than max_burnin")
  counter <- 0
  while(burnin(mp) <= max_burnin & counter < 5){
    message("  k: ", k(hp), ", burnin: ", burnin(mp), ", thin: ", thin(mp))
    mod.list <- replicate(nchains, SingleBatchModel2(dat=dat,
                                                     hp=hp,
                                                     mp=mp))
    mod.list <- harmonizeU(mod.list)
    mod.list <- suppressWarnings(map(mod.list, posteriorSimulation))
    label_swapping <- map_lgl(mod.list, label_switch)
    noswap <- sum(!label_swapping)
    if(noswap < 2){
      burnin(mp) <- as.integer(burnin(mp) * 2)
      mp@thin <- as.integer(thin(mp) * 2)
      ## only increment counter for label switching
      counter <- counter + 1
      mlist <- mcmcList(mod.list)
      neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
      if(is.null(neff))  neff <- 0
      r <- gelman_rubin(mlist, hp)
      next()
    }
    mod.list <- mod.list[ selectModels(mod.list) ]
    mlist <- mcmcList(mod.list)
    neff <- tryCatch(effectiveSize(mlist), error=function(e) NULL)
    if(is.null(neff))  neff <- 0
    r <- gelman_rubin(mlist, hp)
    message("     r: ", round(r$mpsrf, 2))
    message("     eff size (minimum): ", round(min(neff), 1))
    message("     eff size (median): ", round(median(neff), 1))
    if(all(neff > 500) && r$mpsrf < 1.2) break()
    burnin(mp) <- as.integer(burnin(mp) * 2)
    mp@thin <- as.integer(thin(mp) * 2)
    counter <- 0
  }
  model <- combineModels(mod.list)
  meets_conditions <- all(neff > 500) && r$mpsrf < 2 && !label_switch(model)
  if(meets_conditions){
    model <- compute_marginal_lik(model)
  }
  model
})
