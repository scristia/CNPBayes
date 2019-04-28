context("Copy number models")

.test_that(expr, name) NULL

.test_that("Methods defined for the class", {
  sb <- SingleBatchModelExample
  mcmcParams(sb) <- McmcParams(iter=500, burnin=100, thin=1)
  sb <- posteriorSimulation(sb)
  cn.model <- MultiBatchCopyNumber(sb)
  expect_false(manyToOneMapping(cn.model))
  cn <- copyNumber(cn.model)
  zz <- map_z(cn.model)
  expect_identical(as.integer(cn), map_z(cn.model))

  mapping(cn.model) <- rep("1", 3)
  expect_true(manyToOneMapping(cn.model))
  cn.probs <- probCopyNumber(cn.model)
  expect_true(all(as.numeric(cn.probs) == 1))
  cn <- copyNumber(cn.model)
  expect_true(all(cn=="1"))

  mapping(cn.model) <- rep("2", 3)
  cn <- copyNumber(cn.model)
  expect_true(all(cn=="2"))

  mapping(cn.model) <- c("1", "1", "2")
  pz <- probz(cn.model)
  expected <- cbind(rowSums(pz[, 1:2]), pz[, 3])
  cn.probs <- probCopyNumber(cn.model)
  expect_equal(expected, cn.probs)
  cn <- copyNumber(cn.model)
  expected <- z(cn.model)
  expected[expected %in% 1:2] <- "1"
  expected[expected == 3] <- "2"
  expect_identical(cn, expected)

  mapping(cn.model) <- as.character(c(1, 2, 2))
  expected <- z(cn.model)
  expected[expected == 3] <- "2"
  cn <- copyNumber(cn.model)
  expect_identical(cn, expected)
  if(FALSE){
    ## check visualization
    ggMixture(cn.model)
  }
})

.test_that("Mapping components to copy number (single batch)", {
  set.seed(514)
  sb <- SingleBatchModelExample
  cn.model <- SingleBatchCopyNumber(sb)
  params <- mapParams()
  cn.model <- mapComponents(cn.model, params)
  expect_identical(mapping(cn.model), 1:3)
  mapping(cn.model) <- mapping(cn.model)

  map2 <- mapCopyNumber(cn.model, params)
  expect_identical(map2, factor(seq_len(k(cn.model))))

  ##
  ## two components: merge both
  ##
  truth <- simulateData(N=100, p=c(0.9, 0.1),
                        theta=c(0, 0.25), sds=c(0.2, 0.2))
  mp <- McmcParams(iter = 500, burnin = 0, nStarts = 0)
  mcmcParams(truth) <- mp
  model <- posteriorSimulation(truth)

  cn.model <- SingleBatchCopyNumber(model)
  cn.model <- mapComponents(cn.model, params)
  expect_equivalent(mapping(cn.model), c(1L, 1L))
  if(FALSE)
    ggSingleBatch(cn.model)

  cn.map <- mapCopyNumber(cn.model)
  expect_identical(cn.map, factor(rep(2, 2)))
  mapping(cn.model) <- cn.map
  if(FALSE)
    ggSingleBatch(cn.model)
  ##
  ## three components: merge all
  ##
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 0.3), sds=c(0.2, 0.2, 0.2))
  mcmcParams(truth) <- McmcParams(iter=200, burnin=0)
  if(FALSE)
    ggSingleBatch(truth)
  ##model <- posteriorSimulation(truth)
  mp <- McmcParams(iter = 1000, burnin = 100, nStarts = 4)
  hp <- Hyperparameters(tau2.0=100, m2.0=0.1, eta.0=1)
  model <- SingleBatchModel2(mp=mp, dat=y(truth), hp=Hyperparameters(k=3))
  expect_warning(model <- .posteriorSimulation2(model))
  cn.model <- SingleBatchCopyNumber(model)

  cn.model <- mapComponents(cn.model, params)
  expect_equivalent(mapping(cn.model), c(1L, 1L, 1L))
})

.test_that("merge two components", {
  ##
  ## merge 2 of 3 components 
  ##
  set.seed(123)
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 1), sds=c(0.2, 0.2, 0.2))
  mcmcParams(truth) <- McmcParams(iter=200, burnin=0)
  mp <- mcmcParams(truth)
  model <- SingleBatchModel2(mp=mp, dat=y(truth), hp=Hyperparameters(k=3))
  expect_warning(model <- .posteriorSimulation2(model))
  cn.model <- SingleBatchCopyNumber(model)
  ## label switching will occur because components are not well separated
  ## expect_warning(model <- posteriorSimulation(truth),
  ##                "label switching: model k=3")
  cn.model <- SingleBatchCopyNumber(model)
  cn.model <- mapComponents(cn.model, params=mapParams())
  expect_identical(mapping(cn.model), c(1, 1, 2))
  cnmap2 <- mapCopyNumber(cn.model)
  expect_identical(cnmap2, factor(c(2, 2, 3)))
})

.test_that("Mapping components to copy number (multiple batches)", {
  sb <- MultiBatchModelExample
  cn.model <- MultiBatchCopyNumber(sb)
  params <- mapParams()
  cn.model <- mapComponents(cn.model, params)
  expect_identical(mapping(cn.model), 1:3)
})



.test_that("hapmap", {
  set.seed(134)
  dat <- readLocalHapmap()
  b <- collapseBatch(dat, names(dat))
  mp <- McmcParams(iter=1000, burnin=500, nStarts=20)
  ml <- BatchModelList(dat, k=2:5, batch=b, mcmc.params=mp)
  ml <- posteriorSimulation(ml)
  if(FALSE)
    ggMultiBatchChains(ml[[4]])[["batch"]]

  mp <- McmcParams(iter=1000, burnin=500, nStarts=20)
  sb <- SingleBatchModelList(dat, k=2:5, mcmc.params=mp)
  sb <- posteriorSimulation(sb)
  ml.mb <- marginalLikelihood(ml)
  ml.sb <- marginalLikelihood(sb)
  select <- which.max(c(ml.mb, ml.sb))
  ## single-batch 3 components selected
  if(FALSE){
    ggSingleBatch(sb[[2]])
    ggSingleBatch(sb[[4]])
    ggSingleBatch(sb[[5]])
  }



  ggSingleBatch(model)
  ## evaluate merging for k=4
  m4 <- mlist[[3]]
  ggSingleBatch(m4)
  ##
  ## here, component 2 has a large variance
  ##
  ggSingleBatchChains(m4)[["comp"]]





  model <- mlist[[select]]
  d <- densities(model)
  dc <- densitiesCluster(model, merge=TRUE)
  dmlist <- lapply(mlist, DensityModel, merge=TRUE)
  n.comp <- sapply(dmlist, function(x) length(modes(x)))
  ## remove merge models where number components are duplicated
  mlist <- mlist[!duplicated(n.comp)]
  m.y <- marginalLikelihood(mlist)##, params=params)
  argmax <- which.max(m.y)
  expect_true(argmax == 2L)
  if(FALSE){
    plist <- ggSingleBatchChains(mlist[[2]])
    plist[["comp"]]

    plist3 <- ggSingleBatchChains(mlist[[3]])
    plist3[["comp"]]

    ggSingleBatch(mlist[[3]])
    ggSingleBatch(mlist[[2]])

    pstar <- marginal_theta(mlist[[2]])
  }
})
