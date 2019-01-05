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
  library(SummarizedExperiment)
  i <- 1
  data(simulation_parameters, package="PancCnvsData2")
  data(hapmap, package="PancCnvsData2")
  data(cnp_se, package="PancCnvsData2")
  experiment <- simulation_parameters[i, ]
  truth <- hapmap[["truth"]]
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
  mb <- MultiBatchList(data=assays(mb),
                       mp=mcmcParams(mb))
  mlist <- mcmc2(mb)
  expect_identical(modelName(mlist)[1], "MB3")
  p <- ggMixture(mlist[[1]])
  expect_is(p, "gg")
  mb3 <- mlist[[1]]
  expect_identical(mapping(mb3), seq_len(3))
  g <- genotypes(snpdat)[1, ]
  g [ g < 0 ] <- NA
  g <- g+1
  assays(snpdat)[["GT"]][1, ] <- g
  model.lik <- bafLikelihood(mb3, snpdat)
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
  truth2 <- left_join(tibble(id=colnames(snpdat)), truth, by="id")
  ##cnpbayes.stats <- performanceStats(truth2, copyNumber(gt.model))
})

test_that("Compute BAF likelihood for each candidate model", {
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
  mb <- MultiBatch(data=hdat, burnin=150L,
                   iter=200L, max_burnin=150L) %>%
    findSurrogates(0.001) %>% {
    MultiBatchList(data=assays(.),
                   parameters=parameters(.))
    }%>%
    augmentData2
  expect_equal(sum(isSimulated(mb)), 10)
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
  expect_true(validObject(mb2))
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
  mb3 <- mb3[ order(marginal_lik(mb3), decreasing=TRUE) ]
  ## rank by marginal likelihood or proceed?
  expect_identical(k(mb3)[1], 3L)
  tmp <- isAugmented(mb3)
  expect_equal(sum(!tmp), 990L)
  mb4 <- mb3[[1]][!isSimulated(mb3)]
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
  blik <- modelProb(clist[[1]], snpdat2)
  expect_true(blik > -20)
  stats <- baf_loglik(clist, snpdat2)
  expect_identical(stats$cn.model[1], "0,1,2")
  mapping(mb4) <- strsplit(stats$cn.model[1], ",")[[1]]
  model.cn <- factor(copyNumber(mb4))
  pstats <- performanceStats(assays(mb4)$cn, model.cn)
  expect_true(sum(pstats$incorrect) < 3)
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

    library(GGally)
    ggpairs(th)
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

test_that("Rare homozygous deletion exacerbated by batching samples", {
  library(tidyverse)
  set.seed(75)
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
                  c(-1.8, -0.3, 0.1),
                  c(-1.9, -0.4, 0))
  sigmas <- matrix(rep(0.05, 3), 4, 3, byrow=TRUE)
  mb <- simulateBatchData(500, p=ps, theta=thetas, sds=sigmas,
                          batch=dat$batch,
                          zz=dat$z) %>%
    as("MultiBatch")
  assays(mb)$z <- dat$z
  table(batch(mb), assays(mb)$z)
  if(FALSE){
    ggplot(assays(mb), aes(z, oned)) +
      geom_jitter(width=0.1) +
      facet_wrap(~batch)
    ggMixture(mb)
  }
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

