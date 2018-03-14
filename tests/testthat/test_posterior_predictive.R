context("Posterior predictive distribution")

test_that("SB", {
  library(dplyr)
  library(ggplot2)
  set.seed(149)
  tmp <- rst(1000, mean=0, sigma=1)
  expect_equal(mean(tmp), 0, scale=1, tolerance=0.05)
  expect_equal(sd(tmp), 1, scale=1, tolerance=0.05)
  model <- SingleBatchModelExample
  mp <- McmcParams(iter=50, burnin=10)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  tab <- posteriorPredictive(model)
  expect_is(tab, "tbl_df")
  if(FALSE){
    ggPredictive(model, tab)
  }
})

test_that("MB", {
  bmodel <- MultiBatchModelExample
  ##  ch <- chains(bmodel)
  ##  ch <- updateObject(ch)
  ##  ch@u <- matrix(NA, iter(bmodel), k(bmodel))
  ##  chains(bmodel) <- ch
  ##  MultiBatchModelExample <- bmodel
  mp <- McmcParams(iter=10, burnin=5)
  mcmcParams(bmodel) <- mp
  bmodel <- posteriorSimulation(bmodel)
  tab <- posteriorPredictive(bmodel)
  expect_is(tab, "tbl_df")
  expect_identical(colnames(tab), c("y", "component", "batch", "model"))
  if(FALSE){
    ggPredictive(bmodel, tab)
  }
})

test_that("SBP", {
  set.seed(134)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=1)
  sb <- SingleBatchModelExample
  hp <- hpList(k=3)[["SBP"]]
  sbp <- SingleBatchPooled(dat=sample(y(sb)),
                           hp=hp,
                           mp=mp)
  sbp2 <- posteriorSimulation(sbp)
  tab <- posteriorPredictive(sbp2)
  expect_is(tab, "tbl_df")
  if(FALSE) ggPredictive(sbp2, tab)
})

.test_that <- function(nm, expr) NULL

test_that("MBP", {
  set.seed(9)
  mp <- McmcParams(iter=10, burnin=5, nStarts=1, thin=1)
  mb <- MultiBatchModelExample
  mbp.tmp <- as(mb, "MultiBatchPooled")
  expect_identical(length(sigma(mbp.tmp)), 3L)
  ch <- sigma2(chains(mbp.tmp))
  expect_identical(dim(ch), c(iter(mb), 3L))

  hp <- hpList(k=3)[["MBP"]]
  mbp <- MultiBatchPooled(dat=y(mb),
                          hp=hp,
                          mp=mp,
                          batches=batch(mb))
  ch <- sigma2(chains(mbp))
  expect_identical(dim(ch), c(iter(mbp), 3L))
  expect_warning(mbp <- posteriorSimulation(mbp))
  tab <- posteriorPredictive(mbp)
  expect_is(tab, "tbl_df")
  if(FALSE) ggPredictive(mbp, tab)
})

# Plot the posterior predictive distribution for a list of models

.test_that("model lists", {
  library(ggplot2)
  sb <- SingleBatchModelExample
  mcmcParams(sb) <- McmcParams(iter=500, burnin=50)
  sb <- posteriorSimulation(sb)
  models <- list(MultiBatchModelExample, sb)
  tab <- predictiveDataTable(models)
  ggplot(tab, aes(y, fill=predictive)) +
    geom_density(alpha=0.4, adjust=1/2) +
    facet_grid(model~batch) +
    guides(fill=guide_legend(title="")) +
    theme(panel.background=element_rect(fill="white"))
  fig
})
