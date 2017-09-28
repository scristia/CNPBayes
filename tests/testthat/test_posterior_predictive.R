context("Posterior predictive distribution")

test_that("posteriorPredictive", {
  library(dplyr)
  library(ggplot2)
  set.seed(149)
  model <- SingleBatchModelExample
  mp <- McmcParams(iter=500, burnin=50)
  mcmcParams(model) <- mp
  model <- posteriorSimulation(model)
  tab <- posteriorPredictive(model)
  if(FALSE){
    ggPredictive(model, tab)
  }

  bmodel <- MultiBatchModelExample
  mp <- McmcParams(iter=500, burnin=150)
  mcmcParams(bmodel) <- mp
  bmodel <- posteriorSimulation(bmodel)
  tab <- posteriorPredictive(bmodel)
  expect_is(tab, "tbl_df")
  expect_identical(colnames(tab), c("y", "batch", "component"))
  if(FALSE){
    ggPredictive(bmodel, tab)
  }

  if(FALSE){
    set.seed(123)
    ## example from vignette
    se <- readRDS(system.file("extdata", "simulated_se.rds", package="CNPBayes"))
    grl <- readRDS(system.file("extdata", "grl_deletions.rds", package="CNPBayes"))
    cnv.region <- consensusCNP(grl, max.width=5e6)
    i <- subjectHits(findOverlaps(cnv.region, rowRanges(se)))
    med.summary <- matrixStats::colMedians(assays(se)[["cn"]][i, ], na.rm=TRUE)
    mp <- McmcParams(nStarts=5, burnin=1000, iter=1000, thin=2)
    model.list <- gibbs(model="SB", dat=med.summary, k_range=c(2, 3), mp=mp,
                        max_burnin=20000, top=2)
    model <- model.list[[1]]
    tab <- posteriorPredictive(model)
    ggPredictive(model, tab) + xlab("median LRR")

    predict <- tab
    dat <- tibble(y=y(model), predictive="empirical")
    colnames(predict)[2] <- "predictive"

    predict$predictive <- "posterior\npredictive"
    predictive <- NULL
    dat2 <- rbind(dat, predict) %>%
      mutate(predictive=factor(predictive,
                               levels=c("empirical", "posterior\npredictive")))
    fig <- ggplot(dat2, aes(y, fill=predictive)) +
      geom_density(alpha=0.4, adjust=adjust) +
      guides(fill=guide_legend(title="")) +
      theme(panel.background=element_rect(fill="white"))
    fig


    tmp <- filter(dat2, predictive!="empirical")
    table(tmp$y < - 0.2)
    ggplot(tmp,
                  aes(y, ..count.., fill=predictive)) +
      geom_density(alpha=0.4, adjust=adjust) +
      guides(fill=guide_legend(title="")) +
      theme(panel.background=element_rect(fill="white")) +
      xlim(c(-0.6, 0.5))


  }
})
