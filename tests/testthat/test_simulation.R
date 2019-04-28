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

unitTestSimulation <- function(){
  p2 <- 0.01
  twopq <- 2*sqrt(p2)*(1-sqrt(p2))
  q2 <- (1-sqrt(p2))^2
  ps <- matrix(c(p2, twopq, q2),
               4, 3, byrow=TRUE) %>%
    "/"(rowSums(.))
  dat <- tibble(z=sample(1:3, 500,
                         prob=c(p2, twopq, q2), replace=TRUE),
                batch=sort(sample(1:4, 500, prob=rep(1/4, 4),
                                  replace=TRUE)))
  table(dat$z, dat$batch)
  thetas <- rbind(c(-2, -0.5, 0),
                  c(-2.5, -0.7, -0.1),
                  c(-1.8, -0.6, 0.1),
                  c(-1.9, -0.55, 0))
  sigmas <- matrix(rep(0.05, 3), 4, 3, byrow=TRUE)
  mb <- simulateBatchData(500, p=ps,
                          theta=thetas,
                          sds=sigmas,
                          batch=dat$batch,
                          zz=dat$z,
                          df=100) %>%
    as("MultiBatch")
  assays(mb)$z <- dat$z
  mb
}

test_that("Rare homozygous deletion exacerbated by batching samples.", {
  ##
  ## This simulation assumes that we have already identified the batches.  In particular, we use the batches provided by simulateBatchData
  ##
  library(tidyverse)
  set.seed(75)
  mb <- unitTestSimulation()
  truez <- z(mb)
  if(FALSE) ggMixture(mb)
  mb3.all <- MultiBatchList(data=assays(mb))[["MB3"]]
  ##
  ## segfault !
  ##
  burnin(mb3.all) <- 200
  iter(mb3.all) <- 300
  ## this doesn't work
  tmp <- posteriorSimulation(mb3.all)

  ##
  ## Plotting the data, we see homozygous deletion is rare and batch 3 does not have any homozygous deletions
  ##
  ## Excluding all homozygous deletions, we fit the multibatch model with 2 components
  ##
  tmp <- mb[ oned(mb) > -1 ]
  mb2 <- MultiBatchList(data=assays(tmp))[["MB2"]]
  mb2 <- posteriorSimulation(mb2)
  if(FALSE) ggMixture(mb2)
  ##
  ## Fit MB3 for batch 2, using starting values from MB2
  mb3 <- MultiBatchList(data=assays(mb))[["MB3"]]
  sb3 <- mb3[ batch(mb3) == 2 ]
  theta(sb3)[, 2:3] <- theta(mb2)[2, ]
  sigma_(sb3)[, 2:3] <- sigma(mb2)[2, ]
  sigma_(sb3)[, 1] <- rowMeans(sigma_(sb3)[, 2:3, drop=FALSE])
  tau2(sb3) <- rep(tau2(mb2), k(sb3))
  mu(sb3) <- c(NA, mu(mb2))
  p(sb3)[1, ] <- c(sum(oned(sb3) < -1),
              sum(z(mb2)==1),
              sum(z(mb2)==2)) %>%
    "/"(sum(.))
  theta(sb3)[1, 1] <- mean(oned(sb3)[oned(sb3) < -1])
  mu(sb3)[1] <- theta(sb3)[1, 1]
  sigma2.0(sb3) <- sigma2.0(mb2)
  nu.0(sb3) <- nu.0(mb2)
  burnin(sb3) <- 0
  iter(sb3) <- 300
  sb3 <- posteriorSimulation(sb3)
  if(FALSE) ggMixture(sb3)
  ##
  ## Next, simulate homozygous deletions for other batches
  ##
  mb3.all <- MultiBatchList(data=assays(mb))[["MB3"]]
  th <- cbind(theta(sb3)[1, 1], theta(mb2))
  theta(mb3.all) <- th
  s <- cbind(sigma2(sb3)[1, 1], sigma2(mb2))
  sigma2(mb3.all) <- s
  p <- cbind(p(sb3)[1, 1], p(mb2)) %>%
    "/"(rowSums(.))
  p(mb3.all) <- p
  tau2(mb3.all) <- tau2(sb3)
  mu(mb3.all) <- mu(sb3)
  nu.0(mb3.all) <- nu.0(sb3)
  sigma2.0(mb3.all) <- sigma2.0(sb3)
  ## number of batches with missing components
  nMissing <- 3
  ## number observations in missing batches
  nObs <- table(batch(mb3.all))[c(1, 3, 4)] %>%
    as.numeric
  expected_freq <- nObs * p(mb3.all)[c(1, 3, 4), 1]
  expected_freq <- ceiling(expected_freq)

  obsdat <- assays(mb3.all)
  imp <- lapply(expected_freq, function(x, th, s) rnorm(x, th, s),
                th=theta(mb3.all)[c(1, 3, 4), 1],
                s=sigma(mb3.all)[c(1, 3, 4), 1]) %>%
    unlist
  simdat <- tibble(id=paste0("augment_", seq_along(imp)),
                   oned=imp,
                   batch=rep(c(1L, 3L, 4L), expected_freq),
                   is_simulated=TRUE,
                   z=1)
  newdat <- bind_rows(obsdat, simdat) %>%
    arrange(batch)
  ##
  ## Re-initialize model with random starts
  ##
  tmp <- current_values(mb3.all)[c("theta", "sigma2", "mu", "p",
                                   "tau2")]
  mbl <- MultiBatchList(data=newdat,
                        parameters=parameters(mb3.all))
  burnin(100) <- mbl
  iter(0) <- mbl
  posteriorSimulation(mb3.all)

  tmp <- as(mb3.all, "MultiBatchModel")
  logr <- oned(tmp)
  pp <- update_multinomialPr(tmp) %>%
    round(5) %>%
    as.tibble %>%
    mutate(truez=truez,
           logr=round(logr, 3),
           zz=update_z(tmp),
           z2=update_z(tmp))

  p2 <- update_multinomialPr(tmp)
  z2 <- update_z2(p2)
  p2[z2==1, ]
  logr[z2 == 1]
  z3 <- update_z(tmp)

  zz <- update_z(tmp)
  p2 <- update_multinomialPr(tmp)
  z2 <- apply(p2, 1, function(x) sample(1:3, 1, prob=x))

  tmp <- posteriorSimulation(mb3.all)
  ##
  ## 1. Create a MultiBatch object and identify batch surrogates
  ##
  ## this step is already complete.  unitTestSimulation created a MultiBatch object with batch surrogates
  ##

  ##
  ## 2. Assess whether this is a deletion CNP
  ##  - if deletion, run batches with hom. del and batches without hom del. separately.  Or, augment batches without homozygous deletion (?)
  ##
  ##
  object <- mb
  sb <- toSingleBatch(mb)
  iter(sb) <- 200
  burnin(sb) <- 150
  sb <- posteriorSimulation(sb)
  modes(sb) <- computeModes(sb)
  sb <- setModes(sb)
  sds <- sigma(sb)
  fc <- sds/min(sds)
  mn.sd <- c(theta(sb)[1], min(sigma(sb)))

  limits <- mn.sd[1] + c(-1, 1)*2*mn.sd[2]
  freq.del <- assays(object) %>%
    group_by(batch) %>%
    summarize(n = sum(oned < limits[[2]]))
  fewobs <- freq.del$n <= 2
  iszero <- freq.del$n == 0
  expect_false(all(iszero))
  expect_false(!any(fewobs))
  ## at this point, we've established that hom dels are likely in some batches
  ## but not all
  ## -fit k=3 model to batches with homdels
  ## -fit k=2 model to batches with <=2 homdels
  mbl <- MultiBatchList(data=assays(mb))
  iter(mbl) <- 200
  burnin(mbl) <- 150
  mb3 <- mbl[["MB3"]]
  B1 <- filter(freq.del, n >= 1)
  B2 <- filter(freq.del, n < 1)
  mb3 <- mb3[ batch(mb3) %in% B1$batch ]
  mb2 <- mbl[["MB2"]]
  mb2 <- mb2[ batch(mb2) %in% B2$batch ]

  ##dfr(hyn nnperParams(mb3)) <- 10
  mb3 <- posteriorSimulation(mb3)
  ## large variance for one of components in each batch
  ## -- the component with large variance differs between batches
  ## What if there were 10 more homozygous deletions in each of these batches?
  ## ( this is the idea of augmentData2 )
  nsample <- 10
  batches <- rep(unique(batch(mb3)), each=10)
  pred <- predictiveTibble(mb3) %>%
    ##filter(!(batch %in% batches)  & component == 0) %>%
    filter(component == 0) %>%
    "["(sample(seq_len(nrow(.)), sum(nsample), replace=TRUE), ) %>%
    select(oned) %>%
    mutate(batch=rep(batches, nsample),
           id=paste0("augment_", seq_len(nrow(.))),
           oned=oned+rnorm(nrow(.), 0, mn.sd[2])) %>%
    select(c(id, oned, batch))
  newdat <- bind_rows(assays(object), pred) %>%
    arrange(batch) %>%
    mutate(is_simulated = seq_len(nrow(.)) %in% grep("augment_", id))  


  mb2 <- posteriorSimulation(mb2)
  ## - extract thetas and impute missing
  ## - augment data
  ## - fit multibatch model
})

.test_that <- function(expr, name) NULL

.test_that("Simulations in manuscript", {
  ##skip("uses external data")
  library(SummarizedExperiment)
  ##library(pROC)
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
  probe.level <- simulateBatchEffect(hapmap, experiment) %>%
    mutate(is_simulated=FALSE)
  hdat <- meanSummaryBatchEffect(hapmap, probe.level)
  gr <- rowRanges(cnp_se)["CNP_057"]
  snpdat <- hapmapSummarizedExperiment(hapmap, gr)
  snpdat <- snpdat[, hdat$id]
  expect_identical(colnames(snpdat), hdat$id)
  if(FALSE){
    ## sanity check
    ggplot(hdat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  ##
  ## 1. Create a MultiBatch object and identify batch surrogates
  ##
  message("Finding surrogate variables")
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001)

  ##
  ## 2. Assess whether this is a deletion CNP
  ##  - if deletion, run batches with hom. del and batches without hom del. separately.  Or, augment batches without homozygous deletion (?)
  ##

  ##
  ## 3. Create a MultiBatchList of all possible models
  ##
  message("Creating MultiBatchList of all models")
  mb <- mb %>% {
    MultiBatchList(data=assays(.),
                   parameters=parameters(.))
  }
###%>%
  ##augmentData2
  ##expect_equal(sum(isSimulated(mb)), 10)
  snpdat2 <- snpdat[, id(mb)]
  ##snpdat2 <- snpdat[, id2(mb)]
  ##expect_identical(id2(mb), colnames(snpdat2))
  expect_identical(id(mb), colnames(snpdat2))

  ##
  ## 2. Quickly assess which models are most likely to fit this data from an initial burnin
  ##    - sort models by log likelihood
  ##    - use kmeans to cluster the models with the largest log likelihoods
  ##    - select models in the cluster with largest log likelihoods
  message("List the most likely models")
  nStarts(mb) <- 1
  iter(mb) <- 0
  burnin(mb) <- 500
  mb <- lapply(mb, posteriorSimulation)
  loglik <- sapply(mb, log_lik) %>%
    sort(decreasing=TRUE) %>%
    head
  km <- kmeans(loglik, centers=2)
  ix <- order(km$centers, decreasing=TRUE)[1]
  model_names <- names(which(km$cluster == ix))
  mb <- mb[ model_names ]
  mb2 <- listToMultiBatchList(mb)
  ##
  ## 3. For each model, run 2 or more chains starting from the burned in values above
  ## - assess Gelman-Rubin convergence from the multiple chains
  ## - exclude models that have not converged
  nStarts(mb2) <- 2
  iter(mb2) <- 200
  burnin(mb2) <- 100
  expect_true(validObject(mb2))
  expect_identical(theta(mb2[[1]]), theta(mb[[1]]))
  expect_identical(z(mb2[[1]]), z(mb[[1]]))
  for(i in seq_along(mb2)){
    cat(".")
    ## list with length equal to number of starts
    mlist <- as(mb2[[i]], "list")
    m <- posteriorSimulation(mlist) %>%
      combineModels
    ##fails_gr(m) <- is_high_mpsrf(m)
    mb2[[i]] <- m
  }
  not_converged <- !sapply(mb2, convergence)
  mb3 <- mb2[ ! not_converged ]
  ##
  ## 4. For remaining models, compute the marginal likelihood of each
  ##
  iter(mb3) <- 1000
  burnin(mb3) <- 0L
  mb4 <- posteriorSimulation(mb3)
  mb4 <- compute_marginal_lik(mb4)
  ml <- marginal_lik(mb4)
  mb5 <- mb4[ order(ml, decreasing=TRUE) ]
  expect_identical(k(mb5)[1], 3L)
  ##
  ## 5. Select top model
  ##
  topmodel <- mb5[[1]]
  ##  tmp <- isAugmented(mb3)
  ##  expect_equal(sum(!tmp), 990L)
  ##  mb4 <- mb3[[1]][!isSimulated(mb3)]
  ##
  ## 6. Map components to copy number
  ##
  clist <- CnList(topmodel)
  expect_equal(nrow(clist[[1]]), 990)
  if(FALSE){
    ## sanity check
    identical(id(mb5), colnames(snpdat2))
    dat <- tibble(g=genotypes(snpdat2)[1, ],
                  b=bafs(snpdat2)[1, ],
                  cn=copyNumber(clist[[1]]),
                  true.cn=assays(mb5)$cn) %>%
      mutate(cn=factor(cn))
    ggplot(dat, aes(cn, b)) +
      geom_jitter(width=0.1)
    ggplot(dat, aes(true.cn, b)) +
      geom_jitter(width=0.1)
  }
  ##
  ## 7. Calculate probability of the BAFs given the copy number mapping
  ##  - select the mapping that maxizes the probability of the observed allele frequencies
  ##
  blik <- modelProb(clist[[1]], snpdat2)
  expect_true(blik > -20)
  stats <- baf_loglik(clist, snpdat2)
  expect_identical(stats$cn.model[1], "0,1,2")
  mapping(topmodel) <- strsplit(stats$cn.model[1], ",")[[1]]
  ##
  ## 8. For simulation only:  calculate performance characteristics for copy number classification
  ##
  true.cn <- assays(topmodel)$cn
  prob.cn <- probCopyNumber(topmodel)
  cn.score <- prob.cn %*% unique(as.integer(mapping(topmodel)))
  dat <- tibble(true.cn=true.cn,
                score=cn.score[, 1],
                deletion=ifelse(true.cn=="2", "diploid", "deletion"))
  droc <- roc(response=dat$deletion, predictor=dat$score,
              levels=c("deletion", "diploid"),
              auc.polygon=TRUE, plot=FALSE)
  expect_true(droc$auc[[1]] > 0.999)
})


test_that("More challenging simulation", {
  ##skip("uses external data")
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
  experiment <- simulation_parameters[200, ]
  set.seed(132)
  set.seed(experiment$seed[1])
  probe.level <- simulateBatchEffect(hapmap, experiment) %>%
    mutate(is_simulated=FALSE)
  hdat <- meanSummaryBatchEffect(hapmap, probe.level)
  gr <- rowRanges(cnp_se)["CNP_057"]
  snpdat <- hapmapSummarizedExperiment(hapmap, gr)
  snpdat <- snpdat[, hdat$id]
  expect_identical(colnames(snpdat), hdat$id)
  if(FALSE){
    ## sanity check
    ggplot(hdat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001) %>% {
    MultiBatchList(data=assays(.),
                   parameters=parameters(.))
    }%>%
    augmentData2
  snpdat2 <- snpdat[, id2(mb)]
  expect_identical(id2(mb), colnames(snpdat2))
  nStarts(mb) <- 1
  iter(mb) <- 0
  burnin(mb) <- 500
  mb <- lapply(mb, posteriorSimulation)
  loglik <- sapply(mb, log_lik) %>%
    sort(decreasing=TRUE) %>%
    head
  mb2 <- mb[names(loglik)] %>%
    listToMultiBatchList
  nStarts(mb2) <- 3
  iter(mb2) <- 150
  burnin(mb2) <- 300
  for(i in seq_along(mb2)){
    cat(".")
    m <- as(mb2[[i]], "list") %>%
      posteriorSimulation %>%
      combineModels
    fails_gr(m) <- is_high_mpsrf(m)
    mb2[[i]] <- m
  }
  if(FALSE){
    ggMixture(mlist[[1]])
  }
  mb3 <- mb2[ !sapply(mb2, anyWarnings) ] %>%
    compute_marginal_lik
  ## rank by marginal likelihood or proceed?
  mb4 <- mb3[ order(marginal_lik(mb3), decreasing=TRUE) ]
  mb4 <- mb4[[1]][!isSimulated(mb4)]
  clist <- CnList(mb4)
  expect_equal(nrow(clist[[1]]), 990)
  if(FALSE){
    ## sanity check
    identical(id(mb4), colnames(snpdat2))
    dat <- tibble(g=genotypes(snpdat2)[1, ],
                  b=bafs(snpdat2)[1, ],
                  cn=copyNumber(clist[[1]]),
                  true.cn=assays(mb4)$cn) %>%
      mutate(cn=factor(cn))
    ggplot(dat, aes(cn, b)) +
      geom_jitter(width=0.1)
    ggplot(dat, aes(true.cn, b)) +
      geom_jitter(width=0.1)
  }
  stats <- baf_loglik(clist, snpdat2)
  mapping(mb4) <- strsplit(stats$cn.model[1], ",")[[1]]
  model.cn <- factor(copyNumber(mb4))
  ##trace(performanceStats, browser)
  pstats <- performanceStats(assays(mb4)$cn, model.cn)
  expect_true(sum(pstats$incorrect) < 3)
})

test_that("Most challenging simulation", {
  ##skip("uses external data")
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
  experiment <- simulation_parameters[300, ]
  set.seed(132)
  set.seed(experiment$seed[300])
  probe.level <- simulateBatchEffect(hapmap, experiment) %>%
    mutate(is_simulated=FALSE)
  hdat <- meanSummaryBatchEffect(hapmap, probe.level)
  gr <- rowRanges(cnp_se)["CNP_057"]
  snpdat <- hapmapSummarizedExperiment(hapmap, gr)
  snpdat <- snpdat[, hdat$id]
  expect_identical(colnames(snpdat), hdat$id)
  if(FALSE){
    ## sanity check
    ggplot(hdat, aes(cn, baf)) +
      geom_jitter(width=0.1)
  }
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001) %>% {
    MultiBatchList(data=assays(.),
                   parameters=parameters(.))
    }%>%
    augmentTest
  snpdat2 <- snpdat[, id2(mb)]
  expect_identical(id2(mb), colnames(snpdat2))
  nStarts(mb) <- 1
  iter(mb) <- 0
  burnin(mb) <- 500
  mb <- lapply(mb, posteriorSimulation)
  loglik <- sapply(mb, log_lik) %>%
    sort(decreasing=TRUE) %>%
    head
  mb2 <- mb[names(loglik)] %>%
    listToMultiBatchList
  nStarts(mb2) <- 3
  iter(mb2) <- 150
  burnin(mb2) <- 300
  for(i in seq_along(mb2)){
    cat(".")
    m <- as(mb2[[i]], "list") %>%
      posteriorSimulation %>%
      combineModels
    fails_gr(m) <- is_high_mpsrf(m)
    mb2[[i]] <- m
  }
  if(FALSE){
    ggMixture(mlist[[1]])
    tmp <- mb2[["MB3"]]
    tmp <- tmp[!isSimulated(tmp)]
    th <- theta(tmp)
    sds <- sigma(tmp)
    fc <- sds[, 1]/min(sds[, 1])
    th1 <- th[, 1]
    th1[ fc > 2.5 ] <- NA
    th[, 1] <- th1
    yy <- replicate(10, Impute(th), simplify=FALSE)
    th.imputed <- lapply(yy, function(x) x$yimp[, 1]) %>%
      do.call(cbind, .)

    ##library(GGally)
    ##ggpairs(th)
    th2 <- gather(th, "theta", "mean")
    sds2 <- gather(sds, "sigma", "sd")
    library(dplyr)
    mvn <- bind_cols(th2, sds2) %>%
      mutate(component=rep(seq_len(k(tmp)), each=numBatch(tmp)))
  }
  ggMixture(tmp)
  mb3 <- mb2[ !sapply(mb2, anyWarnings) ] %>%
    compute_marginal_lik
  ## rank by marginal likelihood or proceed?
  mb4 <- mb3[ order(marginal_lik(mb3), decreasing=TRUE) ]
  mb4 <- mb4[[1]][!isSimulated(mb4)]
  clist <- CnList(mb4)
  expect_equal(nrow(clist[[1]]), 990)
  if(FALSE){
    ## sanity check
    identical(id(mb4), colnames(snpdat2))
    dat <- tibble(g=genotypes(snpdat2)[1, ],
                  b=bafs(snpdat2)[1, ],
                  cn=copyNumber(clist[[1]]),
                  true.cn=assays(mb4)$cn) %>%
      mutate(cn=factor(cn))
    ggplot(dat, aes(cn, b)) +
      geom_jitter(width=0.1)
    ggplot(dat, aes(true.cn, b)) +
      geom_jitter(width=0.1)
  }
  stats <- baf_loglik(clist, snpdat2)
  mapping(mb4) <- strsplit(stats$cn.model[1], ",")[[1]]
  model.cn <- factor(copyNumber(mb4))
  ##trace(performanceStats, browser)
  pstats <- performanceStats(assays(mb4)$cn, model.cn)
  expect_true(sum(pstats$incorrect) < 3)
})



