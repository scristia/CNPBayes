context("Tile median log R ratios")

test_that("tileMedians", {
  mb <- MultiBatchModelExample
  tiled.medians <- tileMedians(y(mb), 200, batch(mb))
  expected <- tiled.medians %>% group_by(tile) %>%
    summarize(avgLRR=mean(logratio),
              batch=unique(batch))
  tile.summaries <- tileSummaries(tiled.medians)
  expect_identical(tile.summaries, expected)
  mp <- McmcParams(iter=50, burnin=100)
  mb <- MultiBatchModel2(dat=tile.summaries$avgLRR,
                         batches=tile.summaries$batch,
                         mp=mp)
  mb <- posteriorSimulation(mb)
  ggMixture(mb)
  mb2 <- upSample(mb, tiled.medians)
  expect_identical(length(y(mb2)), nrow(tiled.medians))
})
