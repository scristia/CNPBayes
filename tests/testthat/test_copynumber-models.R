context("Copy number models")

  ##
  ## IDEA: There is only one copy number state, but data is not quite normal and
  ## more than a single component is needed to adequately fit the data
  ##
  ##  - the fact that, in truth, there is a single copy number state generating
  ##    should not effect how we fit or select models
  ##
  ##  - copy number inference is a completely separate step that should be
  ##    downstream of model selection
  ##
  ##  - Downstream steps:
  ##       - determine what components belong to a single copy number state
  ##         (a mapping vector of length K denoting distinct copy number states)
  ##           - first step is just to determine which states are distinct
  ##           - second step is to provide the most probable copy number of the
  ##             distinct states
  ##          -  this could be decided based on the degree of overlap or even
  ##             the standard deviation (large variance may indicate outlier component)
  ##       - Extend the SingleBatch and MultiBatch classes to add a @mapping slot
  ##             SingleBatchCopyNumber
  ##             MultiBatchCopyNumber
  ##  - Methods needed for the *BatchCopyNumer classes:
  ##       - k(object)  returns number of distinct copy number states
  ##            - returns number of components if mapping is identity
  ##       - probz(object)
  ##            - when mapping is many to one, posterior probability is added
  ##              for components with the same copy number label
  ##       - gg* plotting methods should be similar, but color code by the mapping
  ##       - Note the marginal likelihood, BIC, etc are defined only for the superclass
  ##           - these methods have nothing to do with the *BatchCopyNumber models
  ##

test_that("Methods defined for the class", {
  sb <- MarginalModelExample
  cn.model <- SingleBatchCopyNumber(sb)

  expect_false(manyToOneMapping(cn.model))
  cn <- copyNumber(cn.model)
  expect_identical(cn, z(cn.model))

  mapping(cn.model) <- rep(1, 3)
  expect_true(manyToOneMapping(cn.model))
  cn.probs <- probCopyNumber(cn.model)
  expect_true(all(as.numeric(cn.probs) == 1))
  cn <- copyNumber(cn.model)
  expect_true(all(cn==1))

  mapping(cn.model) <- rep(2, 3)
  cn <- copyNumber(cn.model)
  expect_true(all(cn==2))

  mapping(cn.model) <- c(1, 1, 2)
  pz <- probz(cn.model)
  expected <- cbind(rowSums(pz[, 1:2]), pz[, 3])
  cn.probs <- probCopyNumber(cn.model)
  expect_equal(expected, cn.probs)
  cn <- copyNumber(cn.model)
  expected <- z(cn.model)
  expected[expected %in% 1:2] <- 1
  expected[expected == 3] <- 2
  expect_identical(cn, expected)

  mapping(cn.model) <- c(1, 2, 2)
  expected <- z(cn.model)
  expected[expected == 3] <- 2
  cn <- copyNumber(cn.model)
  expect_identical(cn, expected)

  if(FALSE){
    ## check visualization
    ggSingleBatch(cn.model)
  }

  mb <- BatchModelExample
  cn.model <- MultiBatchCopyNumber(mb)
  expect_false(manyToOneMapping(cn.model))
  cn <- copyNumber(cn.model)
  expect_identical(cn, z(cn.model))

  mapping(cn.model) <- rep(1, 3)
  expect_true(manyToOneMapping(cn.model))

  cn.probs <- probCopyNumber(cn.model)
  expect_true(all(as.numeric(cn.probs) == 1))
  cn <- copyNumber(cn.model)
  expect_true(all(cn==1))

  mapping(cn.model) <- rep(2, 3)
  cn <- copyNumber(cn.model)
  expect_true(all(cn==2))

  mapping(cn.model) <- c(1, 1, 2)
  pz <- probz(cn.model)
  expected <- cbind(rowSums(pz[, 1:2]), pz[, 3])
  cn.probs <- probCopyNumber(cn.model)
  expect_equal(expected, cn.probs)
  cn <- copyNumber(cn.model)
  expected <- z(cn.model)
  expected[expected %in% 1:2] <- 1
  expected[expected == 3] <- 2
  expect_identical(cn, expected)

  mapping(cn.model) <- c(1, 2, 2)
  expected <- z(cn.model)
  expected[expected == 3] <- 2
  cn <- copyNumber(cn.model)
  expect_identical(cn, expected)

  cn.model <- CNPBayes:::sortComponentLabels(cn.model)
  mapping(cn.model) <- c(1, 2, 2)
  if(FALSE)
    ggMultiBatch(cn.model)
})


test_that("Mapping components to copy number (single batch)", {
  set.seed(514)
  sb <- MarginalModelExample
  cn.model <- SingleBatchCopyNumber(sb)
  params <- mapParams()
  map <- mapComponents(cn.model, params)
  expect_identical(map, 1:3)
  mapping(cn.model) <- map

  map2 <- mapCopyNumber(cn.model, params)
  expect_identical(map2, factor(seq_len(k(cn.model))))

  ##
  ## two components: merge both
  ##
  truth <- simulateData(N=100, p=c(0.9, 0.1),
                        theta=c(0, 0.25), sds=c(0.2, 0.2))
  mp <- McmcParams(iter = 500, burnin = 0, nStarts = 0)
  mcmcParams(truth) <- mp
  expect_warning(model <- posteriorSimulation(truth))

  cn.model <- SingleBatchCopyNumber(model)
  map <- mapComponents(cn.model, params)
  expect_identical(map, c(1L, 1L))
  mapping(cn.model) <- map
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

  model <- posteriorSimulation(truth)
  cn.model <- SingleBatchCopyNumber(model)

  map <- mapComponents(cn.model, params)
  expect_identical(map, c(1L, 1L, 1L))
  mapping(cn.model) <- map

  ##
  ## merge 2 of 3 components 
  ##
  set.seed(123)
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 1), sds=c(0.2, 0.2, 0.2))
  mcmcParams(truth) <- McmcParams(iter=200, burnin=0)
  ## label switching will occur because components are not well separated
  expect_warning(model <- posteriorSimulation(truth),
                 "label switching: model k=3")
  cn.model <- SingleBatchCopyNumber(model)
  map <- mapComponents(cn.model, params)
  mapping(cn.model) <- map
  expect_identical(map, c(1L, 1L, 3L))
  map2 <- mapCopyNumber(cn.model)
  expect_identical(map2, factor(c(2, 2, 3)))
})

test_that("Mapping components to copy number (multiple batches)", {
  sb <- BatchModelExample
  cn.model <- MultiBatchCopyNumber(sb)
  params <- mapParams()
  map <- mapComponents(cn.model, params)
  expect_identical(map, 1:3)

  ##
  ## Assume best fit model is MultiBatch with 3 components. In truth, components
  ## 2 and 3 correspond to 1 copy number state that have more variation than one
  ## would expect if Gaussian.
  ##
  set.seed(100)
  nbatch <- 3
  k <- 3
  means <- matrix(c(-2.1, -1.8, -1.7,
                    -0.4, -0.3, -0.2,
                    -0.45,  -0.25, -0.15), nbatch, k, byrow = FALSE)
  sds <- matrix(0.15, nbatch, k)
  sds[, 1] <- 0.3
  N <- 300
  truth <- simulateBatchData(N = N, batch = rep(letters[1:3],
                                                length.out = N),
                             p = c(1/10, 1/5, 1 - 0.1 - 0.2),
                             theta = means,
                             sds = sds)
  truth@probz[]
  mp <- McmcParams(iter=200, burnin=0, nStarts=0)
  mcmcParams(truth) <- mp
  expect_warning(model <- posteriorSimulation(truth),
                 "label switching: model k=3")
  expect_true(all(range(probz(model)) == c(0, 1)))
  if(FALSE)
    ggMultiBatch(model)

  cn.model <- MultiBatchCopyNumber(model)
  ## trace(CNPBayes:::mapComponents, browser)
  map <- mapComponents(cn.model)
  expect_identical(map, c(1L, 2L, 2L))
  mapping(cn.model) <- map
  ## mapCopyNumber not yet defined for MultiBatch models
  ##cn <- mapCopyNumber(cn.model)
  if(FALSE)
    ggMultiBatch(cn.model)
})



.test_that <- function(expr, name) NULL

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
  sb <- MarginalModelList(dat, k=2:5, mcmc.params=mp)
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

