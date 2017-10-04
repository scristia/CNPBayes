context("Tile median log R ratios")

test_that("tileMedians", {
  mb <- MultiBatchModelExample
  tiled.medians <- tileMedians(y(mb), 200, batch(mb))
  expected <- tiled.medians %>% group_by(tile) %>%
    summarize(avgLRR=mean(logratio),
              batch=unique(batch))
  tile.summaries <- tileSummaries(tiled.medians)
  if(FALSE){
    par(mfrow=c(2, 1))
    hist(y(mb), breaks=100)
    hist(tile.summaries$avgLRR, breaks=100)

    tiled.medians <- tileMedians(y(mb), 20, batch(mb))
    tile.summaries <- tileSummaries(tiled.medians)
    par(mfrow=c(2, 1))
    hist(y(mb), breaks=100)
    hist(tile.summaries$avgLRR, breaks=100)
  }
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
