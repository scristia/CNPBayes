context("MultiBatchList classes")

test_that("Empty MBL", {
  specs <- modelSpecs()
  data <- modelData()
  num_models <- nrow(specs)
  iter <- 1000L; thin <- 1L; nStarts <- 4L
  hp <- Hyperparameters()
  mp <- McmcParams(iter=iter,
                   thin=thin,
                   nStarts=nStarts)
  parameters <- modelParameters(hp=hp, mp=mp)
  chains <- listChains2(specs, parameters)
  current_values <- listValues(specs,
                               data,
                               hp)
  summaries <- listSummaries(specs)
  flags <- listFlags(specs)
  expect_true(validObject(MultiBatchList()))
})

test_that("MBL with data", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  mcmcParams(sb) <- McmcParams(iter=200, thin=2,
                               burnin=200,
                               nStarts=4)
  mb1 <- as(sb, "MultiBatch")
  ## When a single model is not specified,
  ## generate all available models
  mbl <- MultiBatchList(data=assays(mb1))
  expect_true(validObject(mbl))

  ## data with multiple batches
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb1 <- as(mb, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb1))
  L <- length(mbl)
  expect_identical(L, 14L)
  ##
  ## With coercion, we generate a length-1 list
  ##
  mbl <- as(mb1, "MultiBatchList")
  expect_identical(length(mbl), 1L)
  expect_identical(assays(mbl), assays(mb1))
  ##
  ## MultiBatchList -> MultiBatch
  ##  -- only takes first element from list
  mb <- as(mbl, "MultiBatchList")
  expect_true(validObject(mb))

  mbl <- MultiBatchList(data=assays(mb1))
  ## grab a subset of the models
  mbl2 <- mbl[c(1 , 3, 5)]
  expect_identical(specs(mbl2)$model,
                   specs(mbl)$model[c(1, 3, 5)])

  expect_identical(seq_along(mbl), 1:14)
})

test_that("posterior simulation with MBL", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb1 <- as(mb, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb1))
  ## this involves resizing the chains
  mcmcParams(mbl) <- McmcParams(iter=10L,
                                thin=1L,
                                burnin=5L,
                                nStarts=3L)
  expect_true(validObject(mbl))
  ## no code for pooled models yet
  is_pooled <- substr(specs(mbl)$model, 3, 3) == "P"
  mbl <- mbl[ !is_pooled ]
  mbl2 <- posteriorSimulation(mbl)
  expect_true(!any(is.na(theta(chains(mbl2)[[1]]))))
})

test_that("downsampling with lists", {
  ##
  ## with downsampling
  ##
  i <- sort(sample(seq_len(nrow(mbl)), replace=TRUE, size=750))
  ds <- MultiBatchList(data=assays(mb1)[i, ])
  expect_true(validObject(ds))
  iter(ds) <- 100
  burnin(ds) <- 50
  ds <- posteriorSimulation(ds)
  expect_true(validObject(ds))
})

test_that("upsampling with lists", {
  skip("Might be easiest to restrict upsampling to non-list data -- after model selection  ")
  ## Easiest to do upsampling once we have selected a model
  mbl <- full.model
  mbl5 <- upSampleModel(ds, )
  pz <- current_values(mbl5)[[1]][["probz_up"]]
  expect_identical(nrow(pz), specs(mbl5)$number_obs[1])
  expect_false(any(is.na(pz)))
})

## this is a much simpler replacement for gibbs
test_that("mcmc2 for MultiBatch", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  set.seed(1234)
  mb2 <- as(mb, "MultiBatch")
  iter(mb2) <- 100
  burnin(mb2) <- 100
  max_burnin(mb2) <- 200
  nStarts(mb2) <- 4
  ##trace(startingValues2, browser)
  mb.list <- startingValues2(mb2)
  mb3 <- mcmc2(mb2)
  ## select tolerance of 3
  ## - ml should be more stable with more iterations
  expect_equal(-341, marginal_lik(mb3)[[1]], tolerance=3)
})

test_that("augment data for MultiBatch", {
  set.seed(123)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mb <- posteriorSimulation(mb)
  if(FALSE){
    ch.list <- as(chains(mb), "list")
    ggplot(ch.list[["theta"]], aes(s, value, group=b)) +
      geom_line(aes(color=b)) +
      facet_wrap(~k)
  }
  filter <- dplyr::filter
  ## remove homozygous deletions from batch 3
  ## - not a very realistic example because this is not a rare deletion
  dat <- assays(mb) %>%
    filter(! (batch == 3 & oned < -0.5) )
  ## we have problems with this because one batch has only 2 components
  mb1 <- MultiBatch(data=dat)
  iter(mb1) <- 250L
  burnin(mb1) <- 100L
  mb1 <- posteriorSimulation(mb1)
  if(FALSE){
    ## first and second component completely overlap for batch 3
    ggMixture(mb1)
  }
  mbl <- MultiBatchList(data=assays(mb1),
                        mp=mcmcParams(mb1))
  ##
  ## - check for batches with few observations
  ## - augment data if needed
  ## - fit as per usual
  ##
  mbl2 <- augmentData2(mbl)
  expect_true(validObject(mbl2))
  ## 10 observations added
  expect_equal(nrow(mbl2), 1451L, tolerance=3)
  mb3 <- mbl2[[ which(specs(mbl2)$model == "MB3") ]]
  mb3 <- posteriorSimulation(mb3)
  if(FALSE) ggMixture(mb3)
})

test_that("fix probz mcmc2", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  iter(mb) <- 20
  burnin(mb) <- 20
  nStarts(mb) <- 2
  max_burnin(mb) <- 20
  mb2 <- mcmc2(mb)
  expect_true(all(rowSums(probz(mb2))==1))
})

##
## This test is slow
##
test_that("singleBatchGuided", {
  set.seed(1968)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb), burnin=100L,
                        iter=150L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 100L
  sb <- mcmc2(mbl[[ "SB3" ]])
  expect_true( convergence(sb) )

  x <- mbl[ modelName(mbl) != "SB3" ]
  xx <- singleBatchGuided(x, sb)
  expect_identical(names(xx), modelName(x))
  expect_is(xx, "list")
  expect_is(xx[[1]], "list")
  expect_identical(length(xx[[1]]), nStarts(x))
  expect_identical(length(xx), length(x))
  x <- mbl[["MB3"]]
  xx <- singleBatchGuided(x, sb)
  th <- theta(sb)
  expect_equal(colMeans(theta(xx)), th[1, ], tolerance=0.15)
})

test_that("what is wrong with pooled-variance models", {
  set.seed(2000)
  truth <- simulateData(N = 1000,
                        theta = c(-2, -0.4, 0),
                        sds = c(0.05, 0.05, 0.05),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10),
                        df=100)
  mp <- McmcParams(iter = 100, burnin = 100)
  hp <- Hyperparameters(k=3)
  m <- SingleBatchPooled(y(truth), hp, mp)
  model <- posteriorSimulation(m)
  expect_equal(sigma(model)[[1]], 0.05, tolerance=0.02)
  expect_equal(theta(model)[1, ], theta(truth)[1, ], tolerance=0.05)
  if(FALSE) ggMixture(model)

  mb <- as(m, "MultiBatchP")
  mb2 <- posteriorSimulation(mb)
  expect_equal(sigma(mb2)[, 1], 0.05, tolerance=0.02)
  expect_equal(theta(mb2)[1, ], theta(truth)[1, ], tolerance=0.1)
  if(FALSE) ggMixture(mb2)

  ## pooled model with 3 batches
  set.seed(1438)
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  sbp <- as(MultiBatchModelExample, "MultiBatchPooled")
  expect_identical(length(sigma(sbp)), 3L)

  mb <- as(MultiBatchModelExample, "MultiBatchP")
  expect_identical(dim(sigma(mb)), c(3L, 1L))
  nStarts(mb) <- 1L
  mb2 <- posteriorSimulation(mb)
  if(FALSE) ggMixture(mb2)

  ## Try guided function
  mbl <- MultiBatchList(data=assays(mb), burnin=200L,
                        iter=300L)
  mbl <- mbl[ k(mbl) == 3 ]
  max_burnin(mbl) <- 200L
  sb <- mcmc2(mbl[[ "SB3" ]])
  expect_true( convergence(sb) )
  x <- mbl[ modelName(mbl) != "SB3" ]
  ##
  ## Something is wrong with the pooled-variance models
  ## With guided approach, we first fit SB3 (above).
  ## Since SB3 converged, we then fit SBP3, MB3, and MBP3 models
  ## with starting values near the SB3 parameters
  xx <- singleBatchGuided(x, sb)
  sbp <- xx[["SBP3"]][[1]]
  iter(sbp) <- 250L
  burnin(sbp) <- 100L
  sbp <- posteriorSimulation(sbp)
  th <- theta(sbp) %>%
    as.numeric
  expect_equal(th, theta(sb)[1, ], tolerance=0.05)
  if(FALSE) ggMixture(sbp)
  ##
  ## MultiBatch
  ##
  mb <- xx[["MB3"]][[1]]
  mcmcParams(mb) <- mcmcParams(sbp)
  ## This should be faster because we are already in region of high posterior probability
  mb <- posteriorSimulation(mb)
  if(FALSE) ggMixture(mb)
  expect_equal(colMeans(theta(mb)), theta(sb)[1, ],
               tolerance=0.05)
  ##
  ## Guided MCMC2
  ## Now we try the guided version of mcmc2, which basically implements the above steps
  ##mbl <- mbl[ names(mbl) != "SB3" ]
  ##mbl2 <- mcmc2(mbl, sb)
})

test_that("MultiBatchP <-> MultiBatchPooled", {
  x <- MultiBatchPooledExample
  x2 <- as(x, "MultiBatchP") %>%
    as("MultiBatchPooled")
  expect_true(all(names(modes(x)) %in% names(modes(x2))))
})

.test_that <- function(nm, expr) NULL

test_that("Smarter MCMC for MultiBatchList", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ##
  ## beginning of mcmc2
  ##
  object <- MultiBatchList(data=assays(mb), burnin=100L,
                           iter=300L, max_burnin=100L)
  expect_identical(length(unique(batch(object[["SB3"]]))), 1L)
  expect_identical(length(unique(batch(object[["MB3"]]))), 3L)
  sb <- object[["SBP3"]]
  s <- sigma(sb)
  expect_identical(dim(s), c(1L, 1L))
  sb <- object[["MBP3"]]
  expect_is(sb, "MultiBatchP")
  s <- sigma(sb)
  ch <- sigma(chains(sb))
  expect_identical(dim(s), c(3L, 1L))
  expect_identical(dim(summaries(sb)[["data.prec"]]), c(3L, 1L))
  expect_identical(dim(ch), c(iter(sb), nBatch(sb)))

  x <- object[["SBP3"]]
  expect_is(x, "MultiBatchP")
  expect_true(validObject(x))
  x <- object[["MBP3"]]
  expect_is(x, "MultiBatchP")
  expect_true(validObject(x))

  mlist <- listModelsByDecreasingK(object)
  expect_identical(names(mlist), as.character(4:1))
  s <- sigma(mlist[["3"]][["MBP3"]])
  expect_identical(dim(s), c(3L, 1L))
  ##k4 <- fitModelK(mlist[[1]])
  k3 <- fitModelK(mlist[[2]])
  expect_is(k3, "MultiBatchList")
  ## models in list ordered by marginal likelihood
  expect_true(all(convergence(k3)))
  expect_true(all(diff(marginal_lik(k3)) < 0))
  k2 <- fitModelK(mlist[[3]])
  ##expect_false(convergence(k2))

  mlist <- list(k3, k2)
  mlist2 <- as(mlist, "MultiBatchList")
  ix <- order(marginal_lik(mlist2), decreasing=TRUE)
  mlist2 <- mlist2[ix]
  expect_identical(modelName(mlist2)[1], "MB3")

  ##k1 <- fitModelK(mlist[[4]])
  burnin(object) <- 50
  iter(object) <- 100
  max_burnin(object) <- 50
  mlist <- mcmc2(object)
  expect_is(mlist, "MultiBatchList")
  expect_identical(modelName(mlist2)[1], "MB3")
})

test_that("Data not in batch-order", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ix <- sample(seq_len(nrow(mb)), 500, replace=FALSE)
  ##
  ## object not instantiated because data is not in batch order
  ##
  expect_error(MultiBatch(data=assays(mb)[ix, ]))
  expect_error(MultiBatchP(data=assays(mb)[ix, ]))
  expect_error(MultiBatchList(data=assays(mb)[ix, ]))
  mb <- MultiBatchModelExample
  ##
  ## No error is generated with the old classes
  ##
  batch(mb) <- sample(batch(mb), 1500, replace=FALSE)
  validObject(mb)
  expect_error(as(mb, "MultiBatch"))
  expect_error(as(mb, "MultiBatchP"))
  expect_error(as(mb, "MultiBatchList"))
})

.test_that("augment data and down sample for MultiBatchList", {

})

test_that("Generate a list of candidate copy number models", {
  mb <- MultiBatchModelExample
  mb <- as(mb, "MultiBatch")
  mapping(mb) <- seq_len(k(mb))
  expect_identical(mapping(mb), seq_len(k(mb)))
  clist <- CnList(mb)
  expect_is(clist, "CnList")
  mb <- clist[[1]]
  expect_is(mb, "MultiBatch")
  expect_identical(mapping(mb), c("0", "1", "2"))
  expect_identical(cmap(mb), "012")
})

test_that("Compute BAF likelihood for each candidate model", {
  skip("uses external data")
  library(SummarizedExperiment)
  mb <- MultiBatchModelExample
  mb <- as(mb, "MultiBatch")
  mapping(mb) <- seq_len(k(mb))
  clist <- CnList(mb)
  data(cnp_se, package="PancCnvsData2")
  data(hapmap, package="PancCnvsData2")
  data(simulation_parameters, package="PancCnvsData2")
  if(FALSE){
    ## sanity check
    dat <- left_join(hapmap[["b"]], hapmap$truth, by="id")
    ggplot(dat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  experiment <- simulation_parameters[1, ]
  set.seed(132)
  set.seed(experiment$seed[1])
  hdat <- simulateBatchEffect(hapmap, experiment)
  gr <- rowRanges(cnp_se)["CNP_057"]
  snpdat <- hapmapSummarizedExperiment(hapmap, gr)
  snpdat <- snpdat[, hdat$id]
  expect_identical(colnames(snpdat), hdat$id)
  if(FALSE){
    ## sanity check
    ggplot(hdat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  colnames(hdat)[1] <- "provisional_batch"
  colnames(hdat)[6] <- "true_batch"
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001)
  snpdat2 <- snpdat[, id(mb)]
  expect_identical(id(mb), colnames(snpdat2))
  mb <- MultiBatchList(data=assays(mb),
                       mp=mcmcParams(mb),
                       iter=100L,
                       burnin=100L,
                       max_burnin=150L)
  mlist <- mcmc2(mb)
  ix <- grep("augment_", id(mlist))
  expect_identical(assays(mlist), assays(mb))
  if(FALSE){
    ggMixture(mlist[[1]])
  }
  mb3 <- mlist[[1]]
  expect_identical(modelName(mb3), "MB3")
  tmp <- isAugmented(mb3)
  expect_equal(sum(!tmp), 990L)
  clist <- CnList(mb3)
  expect_equal(nrow(clist[[1]]), 990)
  if(FALSE){
    dat <- tibble(g=genotypes(snpdat)[1, ],
                  b=bafs(snpdat)[1, ],
                  cn=copyNumber(clist[[1]]),
                  true.cn=true.z) %>%
      mutate(cn=factor(cn))
    ggplot(dat, aes(cn, b)) +
      geom_jitter(width=0.1)
    ggplot(dat, aes(true.cn, b)) +
      geom_jitter(width=0.1)
  }
  blik <- modelProb(clist[[1]], snpdat2)
  expect_true(blik > -20)
  stats <- baf_loglik(clist, snpdat2)
  ## model fit looks good, but mapping by baf_loglik does not
  expect_identical(stats$cn.model[1], "0,1,2")
  mb3 <- mb3[!isAugmented(mb3), ]
  mapping(mb3) <- strsplit(stats$cn.model[1], ",")[[1]]
  cn <- factor(copyNumber(mb3))
  pstats <- performanceStats(assays(mb3)$cn, cn)
})



