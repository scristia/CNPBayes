context("Copy number models")

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
  sb <- MarginalModelExample
  cn.model <- SingleBatchCopyNumber(sb)
  params <- mapParams()
  map <- mapComponents(cn.model, params)
  expect_identical(map, 1:3)

  ##
  ## two components: merge both
  ##
  truth <- simulateData(N=100, p=c(0.9, 0.1),
                        theta=c(0, 0.25), sds=c(0.2, 0.2))
  mp <- McmcParams(iter = 500, burnin = 0, nStarts = 0)
  mcmcParams(truth) <- mp
  model <- posteriorSimulation(truth)

  cn.model <- SingleBatchCopyNumber(model)
  map <- mapComponents(cn.model, params)
  expect_identical(map, c(1L, 1L))
  mapping(cn.model) <- map
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
  expect_identical(map, c(1L, 1L, 3L))
  mapping(cn.model) <- map
  map <- mapComponents(cn.model, params)
  expect_identical(map, c(1L, 1L, 1L))

  ##
  ## merge 2 of 3 components 
  ##
  truth <- simulateData(N=100, p=c(0.1, 0.8, 0.1),
                        theta=c(-0.3, 0, 1), sds=c(0.2, 0.2, 0.2))
  mcmcParams(truth) <- McmcParams(iter=200, burnin=0)
  model <- posteriorSimulation(truth)
  cn.model <- SingleBatchCopyNumber(model)
  map <- mapComponents(cn.model, params)
  expect_identical(map, c(1L, 1L, 3L))
})

test_that("Mapping components to copy number (multiple batches)", {
  sb <- BatchModelExample
  cn.model <- MultiBatchCopyNumber(sb)
  params <- mapParams()
  map <- mapComponents(cn.model, params)
  expect_identical(map, 1:3)

  ##
  ## Scenario: Suppose best fit model was MultiBatch with 3 components. In
  ## truth, components 2 and 3 correspond to 1 copy number state that have more
  ## variation than one would expect if Gaussian.
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
  mp <- McmcParams(iter=200, burnin=0, nStarts=0)
  mcmcParams(truth) <- mp
  model <- posteriorSimulation(truth)
  if(FALSE)
    ggMultiBatch(model)

  cn.model <- MultiBatchCopyNumber(model)
  map <- mapComponents(cn.model)
  expect_identical(map, c(1L, 2L, 2L))
  mapping(cn.model) <- map
  if(FALSE)
    ggMultiBatch(cn.model)
})
