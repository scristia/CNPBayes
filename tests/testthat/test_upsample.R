context("Down sampling and up sampling")

.test_that <- function(nm, expr) NULL

test_that("upSample2", {
  library(tidyverse)
  set.seed(123)
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1, -0.8, -0.2, 0, 0.2, 0.8, 1, 1.2),
      nbatch, k, byrow = FALSE)
  sds <- matrix(0.1, nbatch, k)
  N <- 1500

  truth <- simulateBatchData(N = N,
                             batch = rep(letters[1:3],
                                         length.out = N),
                             theta = means,
                             sds = sds,
                             p = c(1/5, 1/3, 1 - 1/3 - 1/5))

  ##
  ## Make a tibble:  required plate, plate.index, batch_orig
  ##
  full.data <- tibble(medians=y(truth),
                      plate=batch(truth),
                      batch_orig=as.character(batch(truth))) %>%
    mutate(plate.index=as.integer(factor(plate, levels=unique(plate))))


  ## Below, we down-sample to 500 observations
  ## Required:  batch_orig, batch_index
  partial.data <- downSample(full.data, size=500)
  expect_true(all(c("batch_orig", "batch_index") %in% colnames(partial.data)))

  ##
  ## Required:  a mapping from plate to batch
  ##
  summarize <- dplyr::summarize
  select <- dplyr::select
  plate.mapping <- partial.data %>%
    select(c(plate, batch_index)) %>%
    group_by(plate) %>%
    summarize(batch_index=unique(batch_index))

  ## Fit a model as per usual to the down-sampled data
  mp <- McmcParams(iter=200, burnin=10, )
  hp <- HyperparametersMultiBatch(k=3)
  model <- MultiBatchModel2(dat=partial.data$medians,
                            batches=partial.data$batch_index,
                            mp=mp,
                            hp=hp)
  model <- posteriorSimulation(model)

  ##
  ## Add the batching used for the down-sampled data to the full data
  ##
  full.data2 <- left_join(full.data, plate.mapping, by="plate")
  ##
  ## Estimate probabilities for each individual in the full data
  ##
  model.full <- upSample2(full.data2, model)
  expect_true(validObject(model.full))
  ##
  ## Convert to a CopyNumberModel
  ##
  model.cn <- CopyNumberModel(model.full)
  expect_true(validObject(model.cn))

  ##
  ## An easy way to test whether up-sampling is working, is to not actually
  ## up-sample -- compute the probabilities for the partial data from the
  ## theoretical distributions and compare these probabilties to the posterior
  ## mean of the z's
  ## -- this example is too easy.  Need a more challenging example 
  ##
  model.no.upsampling <- upSample2(partial.data, model, up_sample=FALSE)
  theoretical.p <- probz(model.no.upsampling) %>%
    "/"(rowSums(.))
  model.probs <- probz(model)
  expect_equal(model.probs, theoretical.p, scale=1, tolerance=0.05)


  ## use targeted seq data for a little more challenging
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=500, nStarts=1)
  extfile <- file.path(system.file("extdata", package="CNPBayes"),
                       "targeted_seq.txt")
  dat <- read.delim(extfile)[[1]]
  pdat <- tibble(medians=sample(dat, 500),
                 batch_orig="1",
                 batch="1",
                 batch_index=1L)
  model <- SB(dat=pdat$medians,
              mp=mp,
              hp=hp)
  model <- posteriorSimulation(model)
  model.no.up <- upSample2(pdat, model, up_sample=FALSE)
  theoretical.p <- probz(model.no.up) %>%
    "/"(rowSums(.))
  model.probs <- probz(model)
  expect_equal(as.numeric(theoretical.p),
               as.numeric(model.probs),
               scale=1, tolerance=0.02)
})
