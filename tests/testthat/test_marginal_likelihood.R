context("Marginal likelihood")

.test_that <- function(nm, expr) NULL

.test_that("overfit model", {
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
                        ## default 0.5
                        m2.0=100)
  hp.list <- hpList(mu=-0.75,
                    tau2.0=0.4,
                    eta.0=200, ## default 32
                    ## default 0.5
                    m2.0=100)
  if(FALSE){
    model <- gibbs(model="SB",
                   hp.list=hp.list,
                   mp=mp,
                   dat=galaxies2,
                   batches=rep(1L, length(galaxies2)),
                   k_range=c(2, 3),
                   top=2,
                   max_burnin=8000)
    ##    taut2hat = var(theta(model))
    ## qInverseTau2(mn=0.01, sd=0.001)
    ggChains(model[[1]])
    ggSingleBatch(model[[1]])
    ## The k=3 model has label switching when default hyperparameters are used
    ## But with the hyperparameters in hp.list above, we do not have the label switching issue
    expect_identical(names(model)[1], "SB3")
  }
})

