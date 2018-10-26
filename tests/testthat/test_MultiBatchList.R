context("MultiBatch classes")

test_that("revised_constructors", {
  library(SummarizedExperiment)
  ##
  ##  - get burnin running for multibatch
  ##
  mb <- MultiBatchModel2()
  mcmcParams(mb) <- McmcParams(iter=iter(mb), thin=thin(mb),
                               burnin=burnin(mb))
  expect_true(validObject(mb))
  mp <- McmcParams(iter=10, burnin=1, nStarts=4)
  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample

  ##sb <- updateObject(sb)
  mcmcParams(sb) <- mp
  expect_true(validObject(sb))
  mb2 <- as(sb, "MultiBatch")
  mb2 <- posteriorSimulation(mb2)
  expect_true(validObject(mb2))
  m <- computeModes(mb2)
  dat <- assays(mb2)
  tmp <- MultiBatch("SB3",
                    data=dat,
                    parameters=parameters(mb2))
  identical(parameters(mb2), parameters(tmp))

  ##
  ## By default, we assume no downsampling
  ##
  expect_identical(assays(mb2), downSampledData(mb2))
  ## down sampling is completely random
  down_sample(mb2) <- sample(seq_len(nrow(mb2)), 1000, replace=TRUE)
  expect_error(validObject(mb2))
  ## safer to do downSampleModel
  mb3 <- downSampleModel(mb2, 1000)
  expect_true(validObject(mb3))
  dataMean(mb3)
  ##
  ## MCMC should be run AFTER downsampling.
  ##
  runBurnin(mb3)
  mb3 <- posteriorSimulation(mb3)
  expect_true(validObject(mb3))
  ddat <- downSampledData(mb2)
  expect_identical(nrow(ddat), 1000L)

  ##
  ## Because down_sample for 'mb2' has been specified, the older class automatically receives the downSampledData
  ##
  sb2 <- as(mb3, "MultiBatchModel")
  expect_true(validObject(sb2))
  expect_identical(z(mb3), z(sb2))
  expect_identical(mu(mb3), mu(sb2))
  expect_identical(theta(chains(mb3)), theta(chains(sb2)))
  ##
  ## Turn one specific model (e.g., MB3) into a list with multiple chains
  sb2.list <- as(mb3, "list")
  expect_is(sb2.list[[1]], "MultiBatchModel")
  ## each element of the list should have nStarts set to 1
  expect_identical(nStarts(sb2.list[[1]]), 1L)

  ##
  ## Convert a list with multiple chains back to one single MultiBatch
  ##
  ## old approach
  ##trace(combine_batch, browser)
  tmp <- combine_batch(sb2.list, batches=1L)
  expect_is(tmp, "MultiBatchModel")
  ch <- chains(tmp)
  ch2 <- combineChains(sb2.list)
  expect_identical(ch, ch2)

  combined.mb <- combineModels(sb2.list)
  expect_false(label_switch(combined.mb))
  ## We've combined the 4 chains into a single chain -- nStarts should be 1
  expect_identical(nStarts(combined.mb), 1L)

  ## accessors for values
  theta(mb2)
  sigma2(mb2)
  sigma(mb2)
  p(mb2)
  nu.0(mb2)
  sigma2.0(mb2)
  mu(mb2)
  tau2(mb2)
  logPrior(mb2)
  log_lik(mb2)
  head(probz(mb2))
  head(u(mb2))
  head(z(mb2))

  ## updates for values
  theta(mb2) <- theta(mb2)
  sigma2(mb2) <- sigma2(mb2)
  sigma(mb2) <- sigma(mb2)
  p(mb2) <- p(mb2)
  nu.0(mb2) <- nu.0(mb2)
  sigma2.0(mb2) <- sigma2.0(mb2)
  mu(mb2) <- mu(mb2)
  tau2(mb2) <- tau2(mb2)
  logPrior(mb2) <- logPrior(mb2)
  log_lik(mb2) <- log_lik(mb2)
  probz(mb2) <- probz(mb2)
  u(mb2) <- u(mb2)
  z(mb2) <- z(mb2)

  chains(mb2) <- chains(mb2)
  parameters(mb2) <- parameters(mb2)
  values(mb2) <- values(mb2)
  flags(mb2)
  flags(mb2) <- flags(mb2)

  ## accessors for parameters
  iter(mb2)
  burnin(mb2)
  thin(mb2)
  iter(mb2) <- 10
  burnin(mb2) <- 10
  thin(mb2) <- 10
  hyperParams(mb2)
  mcmcParams(mb2)
  ## replacement
  hyperParams(mb2) <- hyperParams(mb2)
  mcmcParams(mb2) <- mcmcParams(mb2)
  k(mb2)

  ## flags
  label_switch(mb2)
  label_switch(mb2) <- label_switch(mb2)

  ## summaries
  dataMean(mb2)
  dataPrec(mb2)
  dataMean(mb2) <- dataMean(mb2)
  dataPrec(mb2) <- dataPrec(mb2)
  marginal_lik(mb2)
  marginal_lik(mb2) <- marginal_lik(mb2)

  ## data
  assays(mb2)
  oned(mb2)
  batch(mb2)
  oned(mb2) <- oned(mb2)
  ##assays(mb2) <- assays(mb2)
  ##expect_true(validObject(mb2))
})

test_that("starting_values", {
  mb <- MultiBatch()
  expect_true(validObject(mb))

  data(SingleBatchModelExample)
  sb <- SingleBatchModelExample
  dat <- tibble(oned=y(sb),
                batch=1L) %>%
    mutate(id=seq_len(nrow(.)),
           id=as.character(id))

  mb <- MultiBatch(data=dat)
  mb2 <- MultiBatch(data=dat,
                    down_sample=sort(sample(seq_len(nrow(dat)),
                                            500, replace=TRUE)))
  expect_true(validObject(mb2))
})

test_that("findSurrogates", {
  data(SingleBatchModelExample)
  sb <- updateObject(SingleBatchModelExample)
  mb <- as(sb, "MultiBatch")
  assays(mb)[["provisional_batch"]] <- sample(letters[1:15], nrow(mb),
                                              replace=TRUE)
  ##
  ## starting values from modelValues2 are unlikely to be very good
  ##
  mb2 <- findSurrogates(mb)
  expect_true(validObject(mb2))
  iter(mb2) <- 100
  burnin(mb2) <- 50
  thin(mb2) <- 2
  expect_true(validObject(mb2))
  mb2 <- posteriorSimulation(mb2)
})



test_that("downsampling_with_surrogates", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mcmcParams(mb) <- McmcParams(iter=10, thin=2,
                               burnin=50,
                               nStarts=4)
  mb1 <- as(mb, "MultiBatch")
  expect_true(validObject(mb1))
  runBurnin(mb1)
  posteriorSimulation(mb1)

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the chains
  chains(mb1) <- mcmc_chains( specs(mb1), parameters(mb1) )
  mb1 <- posteriorSimulation(mb1)  ## yes

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the summary values
  summaries(mb1) <- modelSummaries( specs(mb1) )
  posteriorSimulation(mb1)  ## yes

  mb1 <- as(mb, "MultiBatch")
  ## will posteriorSimulation work if I reset the current values
  summaries(mb1) <- modelSummaries( specs(mb1) )
  posteriorSimulation(mb1)  ## yes

  mp <- McmcParams(iter=10, burnin=10, thin=1, nStarts=4)
  mcmcParams(mb1) <- mp
  ## data mean / precision should have 3 rows!!
  mb2 <- downSampleModel(mb1, 1000)
  mb2 <- posteriorSimulation(mb2)
  expect_true(validObject(mb2))
  ##
  ## simulate plates from the true batches
  ##
  dat <- assays(mb2)
  dat.list <- split(dat, dat$batch)
  plt.list <- list(letters[1:3], letters[4:6], letters[7:9])
  plt.list2 <- map2(dat.list, plt.list,
                    function(x, plt){
                      sample(plt, nrow(x), replace=TRUE)
                    })
  dat$true_batch <- B <- dat$batch
  pb <- rep(NA, nrow(dat))
  pb[ B==1 ] <- plt.list2[[1]]
  pb[ B==2 ] <- plt.list2[[2]]
  pb[ B==3 ] <- plt.list2[[3]]
  dat$provisional_batch <- pb
  dat %>%
    group_by(provisional_batch) %>%
    summarize(n=n(),
              batch=unique(true_batch))
  assays(mb2) <- dat
  expect_true(validObject(mb2))
  mb3 <- findSurrogates(mb2)
  expect_true(validObject(mb3))
  burnin(mb3) <- 0
  mb4 <- posteriorSimulation(mb3)
})

test_that("upsample", {
  library(SummarizedExperiment)
  library(tidyverse)
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  mcmcParams(mb) <- McmcParams(iter=200, thin=2,
                               burnin=200,
                               nStarts=4)
  mb1 <- as(mb, "MultiBatch")
  mb2 <- downSampleModel(mb1, 1000L)
  burnin(mb2) <- 0L
  mb3 <- posteriorSimulation(mb2)
  expect_identical(specs(mb3)$number_obs,
                   specs(mb1)$number_obs)
  expect_identical(names(values(mb3)),
                   names(modes(mb3)))

  ## downsample with specific index
  i <- sort(sample(seq_len(nrow(mb1)), size=1000L, replace=TRUE))
  mb2 <- downSampleModel(mb1, i=i)
  expect_identical(down_sample(mb2), i)
  ##
  ## with upsampling, we would like to infer the probability of the copy number
  ## assignment from the downsampled data.
  ##
  ## Steps:
  ## 1.  evaluate convergence, etc.
  ## 2.  set current values to modal ordinates
  ## 3.  compute z probabilities for full data from modal ordinates
  ##     (all ordinates are modal except for u)
  ## 4.  give z probabilities from step 3, update z in the standard way
  ## Should we automatically upSample?
  mb4 <- upSampleModel(mb3)
  z1 <- z(mb4)
  z2 <- values(mb4)[["z_up"]]
  z2 <- z2[ down_sample(mb4) ]
  ## check that z is at the modal ordinate for the sampled data
  expect_identical(z1, z2)
})

test_that("Pooled model", {
  mb <- MultiBatchP()
  ## with data
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mb1 <- MultiBatchP(data=assays(mb))
  expect_identical(ncol(sigma(mb1)), 1L)

  tmp <- downSampleModel(mb1)
  tmp2 <- as(tmp, "MultiBatchPooled")

  data(MultiBatchPooledExample)
  mbp <- MultiBatchPooledExample
  mp2 <- as(mbp, "MultiBatchP")
  expect_true(validObject(mp2))
  expect_is(mp2, "MultiBatchP")

  iter(mp2) <- 10L
  burnin(mp2) <- 10L
  mp3 <- posteriorSimulation(mp2)
  log_lik(chains(mp3))
})

test_that("mcmc2 for pooled model", {
  ## with data
  library(SummarizedExperiment); library(tidyverse)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ##
  ## These models should fit the data very well
  ## non-random starts, pooled model
  ## WORKS
  mb1 <- as(mb, "MultiBatchP")
  expect_true(validObject(mb1))
  mb1 <- posteriorSimulation(mb1)
  expect_equal(theta(mb1), theta(mb), tolerance=0.02)
  expect_true(all(rowSums(probz(mb1))== 0))
  ll <- log_lik(chains(mb1))
  expect_identical(length(unique(ll)), iter(mb1))
  mb1 <- posteriorSimulation(mb1)




  ##
  ## random starts, pooled model
  ## ## FAILS
  mb2 <- MultiBatchP(data=assays(mb))
  iter(mb2) <- iter(mb1)
  burnin(mb2) <- burnin(mb1)
  mb2 <- posteriorSimulation(mb2)
  expect_equal(theta(mb2), theta(mb), tolerance=0.02)




  iter(mb1) <- 100L
  burnin(mb1) <- 50L
  nStarts(mb1) <- 2L
  max_burnin(mb1) <- 50L
  mb2 <- mcmc2(mb1)
  expect_is(mb2, "MultiBatchP")

  iter(mb2) <- 500L
  burnin(mb2) <- 1000L
  thin(mb2) <- 2
  mb3 <- posteriorSimulation(mb2)
  if(FALSE) ggMixture(mb3)

  ## this works
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatchPooled")
  mb2 <- posteriorSimulation(mb)
  ggMixture(mb2)
  tmp <- ggChains(mb2)
})

test_that("Plotting", {
  library(tidyverse)
  data(MultiBatchModelExample)
  mbe <- MultiBatchModelExample
  mb <- as(mbe, "MultiBatch")
  ch <- chains(mb)
  ch.list <- as(ch, "list")
  if(FALSE){
    ggplot(ch.list[["theta"]], aes(s, value, group=b)) +
      geom_line(aes(color=b)) +
      facet_wrap(~k)
  }
})

test_that("predictive", {
  library(tidyverse)
  library(magrittr)
  data(MultiBatchModelExample)
  mbe <- MultiBatchModelExample
  prob <- p(mbe)
  z <- seq_len(k(mbe))
  set.seed(123)
  sample(z, replace=TRUE, size=3, prob=prob)
  set.seed(123)
  tmp=replicate(1000, sample_components(z, 3, prob)) %>%
    as.numeric
  phat <- table(tmp)/length(tmp)
  phat <- as.numeric(phat)
  expect_equal(phat, prob, tolerance=0.01)

  set.seed(123)
  update_predictive(mbe)
  ## simulate from rst
  uu <- u(mbe)
  x <- rst(1)

  set.seed(123)
  tmp <- rlocScale_t(1000, 0, 1, 100, uu[1])
  set.seed(123)
  tmp2 <- rst(1000, uu[1], 100, 0, 1)
  expect_identical(tmp, tmp2)

  ch <- chains(mbe)
  predictive(ch) <- matrix(as.numeric(NA),
                           nrow=iter(mbe),
                           ncol=k(mbe)*nBatch(mbe))
  ##
  ## infrastructure set up, but values not quite right
  ##
  ## shouldn't yhat be the same dimension as y?
  ##  i.e., the probability a value is simulated from batch y depends on the representation of the batches.
  chains(mbe) <- ch
  set.seed(123)
  tmp <- cpp_mcmc(mbe, mcmcParams(mbe))
  ystar <- predictive(tmp)
  ##ystar2 <- posteriorPredictive(mbe)
  ##qqplot(as.numeric(ystar), ystar2$y)

  mbe2 <- posteriorSimulation(mbe)
  pmeans <- colMeans(p(chains(mbe2)))
  ystar <- predictive(mbe2)

  mbp <- MultiBatchPooledExample
  mbp2 <- updateObject(mbp)
  expect_true(validObject(mbp2))
  expect_identical(ncol(sigma(chains(mbp2))), nBatch(mbp2))
  expect_identical(ncol(theta(chains(mbp2))), nBatch(mbp2)*k(mbp2))

  tmp <- update_predictiveP(mbp2)
  mbp3 <- posteriorSimulation(mbp2)
  pred <- predictive(mbp3)
  hist(pred, breaks=500)
  hist(y(mbp3), breaks=200)
  p1 <- mean(y(mbp3) < -1)
  p2 <- mean(pred < -1)
  expect_equal(p1, p2, tolerance=0.01)
})


test_that("Empty MBL", {
  specs <- modelSpecs()
  data <- modelData()
  down_sample <- seq_len(nrow(data))
  num_models <- nrow(specs)
  iter <- 1000L; thin <- 1L; nStarts <- 4L
  hp <- Hyperparameters()
  mp <- McmcParams(iter=iter,
                   thin=thin,
                   nStarts=nStarts)
  parameters <- modelParameters(hp=hp, mp=mp)
  chains <- listChains2(specs, parameters)
  current_values <- listValues(specs,
                               data[down_sample, , drop=FALSE],
                               hp)
  summaries <- listSummaries(specs)
  flags <- listFlags(specs)
  MultiBatchList()
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
  expect_identical(L, 15L)
  ##
  ## With coercion, we generate a length-1 list
  ##
  mbl <- as(mb1, "MultiBatchList")
  expect_identical(length(mbl), 1L)
  expect_identical(assays(mbl), assays(mb1))
  expect_identical(down_sample(mbl), down_sample(mb1))
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

  expect_identical(seq_along(mbl), 1:15)
})

## posteriorSimulation now works for MBL models
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
  ##
  ## with downsampling
  ##
  mbl3 <- mbl
  i <- sort(sample(seq_len(nrow(mbl)), replace=TRUE, size=750))
  tmp <- downSampleModel(mbl[[5]], i=i)
  expect_true(validObject(tmp))

  mbl3 <- downSampleModel(mbl, N=750)
  expect_identical(length(u(mbl3[[1]])), 750L)
  mb <- mbl3[[1]]
  expect_true(validObject(mb))
  posteriorSimulation(mb)
  expect_true(all(is.na(theta(chains(mbl3)[[1]]))))
  mbl4 <- posteriorSimulation(mbl3)
  expect_true(all(!is.na(theta(chains(mbl4)[[1]]))))

  mbl5 <- upSampleModel(mbl4)
  pz <- values(mbl5)[[1]][["probz_up"]]
  expect_identical(nrow(pz), specs(mbl5)$number_obs[1])
  expect_false(any(is.na(pz)))
})

## this is a much simpler replacement for gibbs
test_that("mcmc2 for MultiBatch", {
  data(MultiBatchModelExample)
  mb <- MultiBatchModelExample
  set.seed(123)
  mb2 <- as(mb, "MultiBatch")
  mb3 <- mcmc2(mb2)
  expect_equal(-341, marginal_lik(mb3)[[1]], tolerance=1)
})

test_that("mcmc2 for MultiBatchList", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  mbl <- MultiBatchList(data=assays(mb))
  mp <- McmcParams(burnin=50, iter=100, nStarts=4,
                   max_burnin=100)
  mcmcParams(mbl) <- mp
  set.seed(123)
  mbl2 <- mcmc2(mbl)
  ml <- marginal_lik(mbl2)
  expect_true(all(diff(ml) <= 0))
})

test_that("augment data for MultiBatch", {
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
  mb1 <- posteriorSimulation(mb1)
  ## this looks terrible!
  if(FALSE) ggMixture(mb1)
  mbl <- MultiBatchList(data=assays(mb1))
  burnin(mbl) <- 200

  ##
  ## - check for batches with few observations
  ## - augment data if needed
  ## - fit as per usual
  ##
  mbl2 <- augmentData2(mbl)
  expect_true(validObject(mbl2))
  ## 10 observations added
  expect_identical(nrow(mbl2), 1454L)
  mb3 <- mbl2[[ which(specs(mbl2)$model == "MB3") ]]
  expect_identical(seq_len(nrow(mb3)), down_sample(mb3))
  mb3 <- posteriorSimulation(mb3)
  if(FALSE) ggMixture(mb3)
})

test_that("fix probz in mcmc2", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  iter(mb) <- 100
  burnin(mb) <- 50
  nStarts(mb) <- 4
  mb2 <- mcmc2(mb)
  expect_true(!any(is.na(probz(mb2))))

  ## try mcmc2 for multibatch list
  mbl <- MultiBatchList(data=assays(mb2))
  mbl <- mbl[1:4]
  mcmcParams(mbl) <- mcmcParams(mb2)
  set.seed(123)
  mbl2 <- mcmc2(mbl)
  expect_true(all(rowSums(probz(mbl2[[1]])) == 1))
})

test_that("Smarter MCMC for MultiBatchList", {
  library(SummarizedExperiment)
  data(MultiBatchModelExample)
  mb <- as(MultiBatchModelExample, "MultiBatch")
  ##
  ## beginning of mcmc2
  ##
  object <- MultiBatchList(data=assays(mb), burnin=250L)
  x <- object[["SBP3"]]
  expect_is(x, "MultiBatchP")
  expect_true(validObject(x))
  expect_identical(250L, burnin(object))
  N <- nrow(object)
  object2 <- augmentData2(object)
  mp <- mcmcParams(object2)
  max_burnin(mp) <- 100
  iter(mp) <- 250
  mcmcParams(object2) <- mp
  ## no data should be added -- this is a common deletion
  expect_identical(N, nrow(object2))
  sp <- specs(object2)
  object2.list <- split(object2, sp$k)
  ##
  ## reverse the order so that model with largest number of components is first
  ##
  object2.list <- rev(object2.list)
  ##
  ## For models with k > 1, fit the SB model first.
  ##
  ## If label switching occurs, do not fit the other models.
  ##
  ##
  model.nms <- sapply(object2.list, k) %>%
    sapply("[", 1) %>%
    paste0("SB", .)
  for(i in seq_along(object2.list)){
    model.list <- object2.list[[i]]
    sb <- model.list[[1]]
    mod.list <- model.list[-1]
    ##
    ## only fit multibatch models for given k if the corresponding single-batch model converges
    ##
    sb2 <- mcmc2( sb )
    if( convergence(sb2) ){
      for(j in seq_along(mod.list)){
        ##singleBatchGuided(mod.list[[j]], sb2)
        x <- mod.list[[j]]
        expect_is(x, "MultiBatchP")
      }
    } else {
      ## only keep the single-batch model
      model.list <- model.list[1]
    }
    object2.list[[i]] <- model.list
  }


    object2.list[[i]][[sb.model]] <- SB2
    if( convergence(SB2) ){
      ## evaluate other models in this list
      ## simulate starting values from modal values

    } else{
      object2.list[[i]] <- object2.list[[i]][ sb.model ]
    }
  }




    pred2 <- pred %>%
      longFormatKB(K=k(object), B=numBatch(object)) %>%
      set_colnames(c("s", "oned", "batch", "component")) %>%
      mutate(model=modelName(object)) %>%
      mutate(batch=as.character(batch),
             batch=gsub("batch ", "", batch)) %>%
      select(-component) ## component labels are wrong
    zz <- zstar(object) %>%
      longFormatKB(K=k(object), B=numBatch(object)) %>%
      mutate(component=factor(value))
    pred2$component <- zz$component
    ##
    ## if not converged, skip remaining models in this category
    ##
##  }

  ix <- order(sp$k, decreasing=TRUE)
  object2 <- object2[ix]

  for(i in seq_along(object)){
    object[[i]] <- mcmc2(object[[i]])
  }
  ml <- marginal_lik(object)
  object <- object[ !is.na(ml) ]
  ml <- ml [ !is.na(ml) ]
  object <- object[ order(ml, decreasing=TRUE) ]
  object
})

test_that("augment data and down sample for MultiBatchList", {

})

test_that("genotype mixture components MultiBatch", {

})

test_that("genotype mixture components MultiBatchList", {

test_that("mcmc2 with MBL", {
  

})
