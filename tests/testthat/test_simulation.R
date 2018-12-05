context("Simulation")

test_that("test_simulation", {
  set.seed(42)
  arguments <- list(sl.good = 6.25, sl.bad = 0.0625, prbias = 0.03,
                    n = 0.2, prvar = c(19.92985, 0.06272))
  dat <- simulateProbeLevel(cnvs = 1, K = 4, probes = 10,
                            arguments = arguments, qual = "medium")
  expect_identical(names(dat), c("measurements", "assignments"))
  if(FALSE){
    x <- dat[[1]]
    if (FALSE)
      hist(rowMeans(x[, , 1, 3]), col = "gray", breaks = 80)
    K <- 4
    xx <- x[, , 1, K]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center = TRUE, scale. = TRUE)$x[, 1]
    if (cor(pc, mns) < cor(-pc, mns)) {
      pc <- -pc
    }
    mp <- McmcParams(iter = 100, burnin = 10, nStarts = 1)
    mlist <- SingleBatchModelList(data=pc, k=1:4, mcmc.params=mp)
    mlist <- posteriorSimulation(mlist)
    m.y <- marginalLikelihood(mlist)
    expect_true(which.max(m.y) == 4L)
    if (FALSE){
      hist(pc, breaks = 100, col = "gray", border = "gray")
    }
  }
})

test_that("simulations in manuscript", {
  skip("test uses external data")
  ##load_all("../../../PancCnvs")
  library(PancCnvsData2)
  library(tidyverse)

  i <- 1
  hapmap.dir <- system.file("extdata/hapmap", package="PancCnvsData2")
  cnp_57 <- readRDS(file.path(hapmap.dir, "cnp_57.rds"))
  true.z <- readRDS(file.path(hapmap.dir, "cnp57_truth.rds"))

  expect_identical(names(cnp_57), c("r", "b", "g"))
  hdat <- hapmapData(cnp_57, true.z, experiment[i, ])
  set.seed(experiment$seed[i])
  mb <- MultiBatch(data=hdat, burnin=100L,
                   iter=300L, max_burnin=100L)
  mb2 <- findSurrogates(mb, 0.001)
  mb <- MultiBatchList(data=assays(mb2),
                       mp=mcmcParams(mb2))
  mlist <- mcmc2(mb)
  expect_identical(modelName(mlist)[1], "MB3")
  p <- ggMixture(mlist[[1]])
  mb3 <- mlist[[1]]
  expect_identical(mapping(mb3), seq_len(3))

  model.lik <- bafLikelihood(mb3, snpdat)

  expect_is(p, "gg")
  if(FALSE){
    mp <- McmcParams(iter=1000, burnin=5000, nStarts=4, thin=2)
    models <- gibbs(model=c("SB", "MB", "MBP", "SBP"),
                    dat=full.data$medians,
                    batch=full.data$batch_index,
                    k_range=c(1, 4),
                    mp=mp,
                    max_burnin=5000,
                    df=100,
                    min_effsize=200)
  }

  snpdat <- hapmapSummarizedExperiment(cnp_57)
  names(true.z) <- colnames(snpdat)
  dnms <- list(NULL, c("sensitivity", "specificity"))
  cn.model <- CopyNumberModel(mlist[[1]])

  if(k(cn.model) == 1){
    mapping(cn.model) <- "2"
    gt.model <- cn.model
    loglik <- tibble(model="2", loglik=0)
  } else {

    gt.model <- model.lik[["model"]]
    loglik <- model.lik$loglik
  }
  truth2 <- truth[colnames(snpdat)]
  cnpbayes.stats <- performanceStats(truth2, copyNumber(gt.model))
})

test_that("more challenging simulation", {
  skip("more challenging simulation")
  library(SummarizedExperiment)
  library(magrittr)
  library(devtools)
  library(tidyverse)
  i <- 200
  data(simulation_parameters, package="PancCnvsData2")
  data(hapmap, package="PancCnvsData2")
  data(cnp_se, package="PancCnvsData2")
  experiment <- simulation_parameters[i, ]
  set.seed(experiment$seed)
  fname <- paste0("shift_", i, ".rds")
  hdat <- simulateBatchEffect(hapmap, experiment)
  snpdat <- hapmapSummarizedExperiment(hapmap,
                                       rowRanges(cnp_se)["CNP_057"])
  colnames(hdat)[1] <- "provisional_batch"
  colnames(hdat)[6] <- "true_batch"
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.0001)
  ## batch 4 has different mixture proportions
  mbl <- MultiBatchList(data=assays(mb),
                        iter=100L,
                        burnin=500L,
                        max_burnin=150L,
                        nStarts=4)
  mb <- mbl[["MB3"]]
  mb <- mcmc2(mb)
  ##summaries(mb)[["data.mean"]] <- computeMeans(mb)
  ##summaries(mb)[["data.prec"]] <- computePrec(mb)
  mb <- mcmc2(mbl[["MB3"]])



  tmp <- posteriorSimulation(mb)
  mb4 <- mb[batch(mb) == 4]
  mb4 <- posteriorSimulation(mb4)
  tmp2 <- mb[ batch(mb) != 4]
  tmp2 <- posteriorSimulation(tmp2)
  tmp2 <- mb[ !batch(mb) %in% c(3, 4) ]
  tmp2 <- posteriorSimulation(tmp2)
  counts=tableBatchZ(mb2)
  batch_totals=rowSums(counts)
  weighted_sums <- rep(NA, k(mb2))
  for(k in 1:3){
    weighted_sums[k] <- sum(counts[, k] * batch_totals) / sum(batch_totals)
  }

  ## problem batch

  ## TODO: I don't like this, but we currently require batches to be renumbered from 1
  ## So batch 4 is now batch 1
  expect_true(all(batch(mb4) == 1))
  mb4 <- posteriorSimulation(mb4)
  ## fit for this batch is perfectly reasonable
  ggMixture(mb4)

  sigma(tmp)
  theta(tmp)[4, ] <- colMeans(theta(tmp))
  sigma2(tmp)[4, ] <- colMedians(sigma2(tmp))

  expect_true(validObject(mb4))
  mb4 <- posteriorSimulation(mb4)

  snpdat2 <- snpdat[, id(mb)]
  mb <- MultiBatchList(data=assays(mb),
                       mp=mcmcParams(mb),
                       iter=100L,
                       burnin=100L,
                       max_burnin=150L)
  ## short MCMC with large number of candidate models
  mlist <- mcmc2(mb)
  nStarts(mlist) <- 4
  iter(mlist) <- 1000
  burnin(mlist) <- 500
  iter(mlist) <- 300
  burnin(mlist) <- 300
  mlist2 <- mcmc3(mlist)
  mb <- mlist2[[1]]
  ix <- grep("augment_", id(mb))
  clist <- CnList(mb)
  stats <- baf_loglik(clist, snpdat2)
  mapping(mb) <- strsplit(stats$cn.model[1], ",")[[1]]
  cn <- factor(copyNumber(mb))
  pstats <- performanceStats(as.integer(assays(mb)$cn), as.integer(cn))
}
