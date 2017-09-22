context("Down sampling")

## For large datasets, downsample observations. This makes it easier to run a
## large number of MCMC simulations for big datasets.

.test_that <- function(name, expr) NULL

.test_that("downsample", {
  ##meds <- colMedians(xx)
  dat <- readLocalHapmap()
  mp <- McmcParams(iter=1000, burnin=1000,
                   nStarts=4, thin=10)
  batches <- collapseBatch(dat, names(dat), THR=0.05) %>%
    factor %>% as.integer
  hp.list <- hpList(k=5)
  set.seed(123)
  sb <- SingleBatchPooled(dat=dat, hp=hp.list[["SBP"]], mp=mp)
  sb <- posteriorSimulation(sb)
  ggMixture(sb)

  ##
  ## downsample
  ##
  ds <- downSampleEachBatch(dat, 200, batch=batches)
  sb2 <- SingleBatchPooled(dat=ds$y, hp=hp.list[["SBP"]], mp=mp)
  sb2 <- posteriorSimulation(sb2)
  ggMixture(sb2)

  sb3 <- upSample(sb2, ds)
  expect_identical(sort(y(sb)), sort(y(sb3)))
  pz <- probz(sb3)
  expect_equal(range(pz), c(0, 1))
  expect_identical(nrow(pz), length(y(sb)))
  if(FALSE){
    ggMixture(sb3)
  }
  ##
  ## pseudocode
  ##
  if(FALSE){
    model.fulldata <- SingleBatchModel(data)
    model.partdata <- downSample(model.fulldata)
    model.partdata <- posteriorSimulation(model.partdata)
    model.fulldata1 <- upSample(model.partdata)
    model.fulldata2 <- posteriorSimulation(model.fulldata)
    ## compare z's
    ## compare thetas
    ## compare likelihood
    identical(model.fulldata1, model.fulldata2)
  }
})
