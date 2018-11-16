context("Down sampling")
test_that("down sampling", {
  library(tidyverse)
  set.seed(1)
  truth <- simulateData(N = 5000,
                        p = rep(1/3, 3),
                        theta = c(-1, 0, 1),
                        sds = rep(0.1, 3),
                        df=100)
  mp <- McmcParams(iter = 100, burnin = 100)
  model <- SingleBatchModel2(dat = y(truth),
                             hp=hpList(k=3)[["SB"]],
                             mp = mp)
  ## simulate 50 provisional batch labels
  plates <- expand.grid(letters[1:26], letters[1:26]) %>%
    apply(1, paste0, collapse="") %>%
    unique %>%
    "["(1:50)
  plates <- sample(plates, 5000, replace=TRUE)
  plate.index <- as.numeric(factor(plates))
  full.data <- tibble(medians=y(truth),
                      provisional_batch=plates) %>%
    mutate(provisional.index=plate.index) %>%
    mutate(batch=collapseBatch(medians, plate.index, 0.05))
  full.data2 <- full.data %>%
    mutate(batch_orig=factor(as.integer(factor(batch)))) %>%
    mutate(batch_orig=as.character(batch_orig))
  ## Downsample to 750 observations
  partial.data <- downSample(full.data2,
                             size=750,
                             min.batchsize=50)
  ## check that all batches were sampled
  iter <- 0
  L1 <- length(unique(partial.data$provisional_batch))
  L2 <- length(unique(full.data2$provisional_batch))
  size <- 750
  while(L1 < L2 && iter < 3){
    partial.data <- downSample(full.data2, size=size)
    L1 <- length(unique(partial.data$provisional_batch))
    iter <- iter+1
  }
  partial.data
  expect_identical(nrow(partial.data), 750L)
  ## since no batch effect was simulated, expect few batches in the downsampled data
  expect_true(length(unique(partial.data$batch)) < 3)
})
