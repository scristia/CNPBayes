context("analysis of rare deletions")

## for first step, see test_summarized_region.R
.test_that <- function(nm, expr) NULL

test_that("augment data", {
  set.seed(123)
  seeds <- sample(seq_len(10000), nrow(cnp_se), replace=TRUE)
  set.seed(seeds[ 1 ])
  THR <- -1
  mb.subsamp <- readRDS("../../inst/extdata/mb_subsamp.rds")
  ##
  ## With a bad choice of start, we can get stuck in a local mode for a very long time. Look at several starts and select the most promising among these.
  ##
  ## - 10 random starts for both SB3 and SBP3 models
  ## - additional MCMC simulations for model with highest log lik
  ##
  if(sum(oned(mb.subsamp) < THR) < 5){
    message("Augment data for rare deletions (less than 5 obs)")
    hdmean <- mean(oned(mb.subsamp)[ oned(mb.subsamp) < THR])
    hdsd <- sd(oned(mb.subsamp)[ oned(mb.subsamp) < THR])
    if(is.na(hdsd)) hdsd <- 0.1
    loc.scale <- tibble(theta=hdmean,
                        sigma2=hdsd^2,
                        phat=0.01,
                        batch=seq_len(nrow(theta(mb.subsamp))))
    imp.hd <- impute(mb.subsamp, loc.scale, start.index=1)
    obsdat <- assays(mb.subsamp) %>%
      mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat, imp.hd) %>%
      arrange(batch)
  } else {
    simdat <- assays(mb.subsamp)
  }
  if(FALSE){
    saveRDS(simdat, file="../../inst/extdata/simdat.rds")
  }
  expected <- readRDS("../../inst/extdata/simdat.rds")
  expect_equivalent(simdat, expected)
})

test_that("warmup", {
  set.seed(5)
  simdat <- readRDS("../../inst/extdata/simdat.rds")
  ##
  message("Random starts and short warmup for SBP3 model")
  ##
  mbl <- replicate(10, MultiBatchList(data=simdat)[["SBP3"]])
  for(j in seq_along(mbl)){
    cat(".")
    mb <- mbl[[j]]
    iter(mb) <- 0
    burnin(mb) <- 100
    mb <- posteriorSimulation(mb)
    mbl[[j]] <- mb
  }
  ##
  message("Random starts and short warmup for SB3 model")
  ##
  ml <- sapply(mbl, log_lik)
  mbl2 <- replicate(10, MultiBatchList(data=simdat)[["SB3"]])
  for(j in seq_along(mbl)){
    cat(".")
    mb <- mbl2[[j]]
    iter(mb) <- 0
    burnin(mb) <- 100
    mb <- posteriorSimulation(mb)
    mbl2[[j]] <- mb
  }
  ##
  message("Choose SB3 or SBP3 model with highest log likelihood.  Penalize SB3 model for the additional parameters by adding 50.")
  ##
  ml2 <- sapply(mbl2, log_lik)
  if(max(ml2, na.rm=TRUE) > max(ml, na.rm=TRUE) + 50){
    sb3 <- mbl2[[which.max(ml2)]]
  } else {
    sb3 <- mbl[[which.max(ml)]]
  }
  if(FALSE){
    saveRDS(sb3, file=file.path("..",
                                "..",
                                "inst",
                                "extdata",
                                "sb3.rds"))
  }
  expected <- readRDS(file.path("..",
                                "..",
                                "inst",
                                "extdata",
                                "sb3.rds"))
  expect_equivalent(sb3, expected)
})

test_that("mcmc1", {
  set.seed(123)
  seeds <- sample(seq_len(10000), nrow(cnp_se), replace=TRUE)
  set.seed(seeds[ 1 ])
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3.rds"))
  ##
  message("Run 500 burnin and 400 additional simulations.")
  ##
  iter(sb3) <- 400
  burnin(sb3) <- 500
  sb3 <- posteriorSimulation(sb3)
  expect_equal(as.numeric(theta(sb3)),
               c(-3.76, -0.217, -0.0019),
               tolerance=0.01)
  if(FALSE){
    saveRDS(sb3, file=file.path("..",
                                "..",
                                "inst",
                                "extdata",
                                "sb3_1.rds"))
  }
  expected <- readRDS(file.path("..",
                                "..",
                                "inst",
                                "extdata",
                                "sb3_1.rds"))
  expect_equivalent(sb3, expected)
})

test_that("early stop", {
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3_1.rds"))
  pz <- probz(sb3)
  maxprob <- rowMax(pz)
  (number_lowprob <- sum(maxprob < 0.95))
  pz <- probz(sb3) %>%
    "/"(rowSums(.)) %>%
    rowMax
  ##
  message("Compute fraction of samples with posterior probability > 0.99")
  mean_maxp <- mean(pz > 0.99)
  expect_equal(mean_maxp, 0.888, tolerance=0.001)
  ##nif(mean_maxp > 0.98 || CNP == "CNP_058"){
  condition <- mean_maxp > 0.995
  expect_false(condition)
})

test_that("Fit multiple batches.  First, pooled-variance model without homozygous deletions", {
  ##
  ## Read back in the data with all batches
  ##
  mb.subsamp <- readRDS(file.path("..",
                                  "..",
                                  "inst",
                                  "extdata",
                                  "mb_subsamp.rds"))

  THR <- -1
  ##
  ##
  ##
  message("Fit MBP2 model excluding homozygous deletions")
  ##
  ##
  ##
  fdat <- filter(assays(mb.subsamp), oned > THR)
  mbl <- replicate(10, MultiBatchList(data=fdat)[["MBP2"]])
  for(j in seq_along(mbl)){
    cat(".")
    mb <- mbl[[j]]
    burnin(mb) <- 100
    iter(mb) <- 0
    mb <- posteriorSimulation(mb)
    mbl[[j]] <- mb
  }
  mb <- mbl[[ which.max(sapply(mbl, log_lik)) ]]
  iter(mb) <- 500
  burnin(mb) <- 50
  mod_2.3 <- posteriorSimulation(mb)
  assays(mod_2.3)$is_simulated=FALSE
  if(FALSE){
    saveRDS(mod_2.3, file=file.path("..",
                                    "..",
                                    "inst",
                                    "extdata",
                                    "mod_2.3.rds"))
  }
  expected <- readRDS(file.path("..",
                                "..",
                                "inst",
                                "extdata",
                                "mod_2.3.rds"))
  expect_equivalent(mod_2.3, expected)
})

test_that("Compare variances in pooled variance model to variance from SingleBatch model", {
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3_1.rds"))
  varratio <- max(sigma2(sb3))/min(sigma2(sb3))
  expect_equal(varratio, 84.38, tolerance=0.01)
  augment_hemizygous <- p(sb3)[2] < 0.1 || varratio > 100
  expect_true(augment_hemizygous)
})

test_that("Augment hemizygous", {
  ##
  ## Reminder that we only do this section if augment_hemizygous
  ## is TRUE
  ##
  mod_2.3 <- readRDS(file.path("..",
                               "..",
                               "inst",
                               "extdata",
                               "mod_2.3.rds"))
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3_1.rds"))
  augment_hemizygous <- TRUE
  if(!augment_hemizygous) stop()
  is_pooledvar <- TRUE
  batch_labels <- assays(mod_2.3)$batch_labels
  batch_labels <- batch_labels %>%
    factor(., levels=unique(.)) %>%
    as.integer(.) %>%
    unique(.) %>%
    sort()
  ##
  message("Data augmentation for hemizygous deletion component")
  ##
  ##
  ## normalize probabilities
  pz <- probz(mod_2.3) %>%
    "/"(rowSums(.)) %>%
    rowMax
  ##
  ## batches with high posterior probabilities
  ##
  ubatch <- unique(batch(mod_2.3)[ pz >= 0.95 ])
  ##
  ## Find number of samples assigned to the hemizygous deletion
  ## component with high probability.
  ##
  ## Check if any of the batches have fewer than 10 subjects
  ## with high posterior probability
  ##
  tmp <- tibble(batch=batch(mod_2.3),
                z=map_z(mod_2.3),
                pz=pz) %>%
    group_by(batch) %>%
    summarize(N=n(),
              n=sum(z==1 & pz > 0.9)) %>%
    filter(n < 10)
  ##
  ##
  ##
  is_dropped <- !batch_labels %in% ubatch |
    batch_labels %in% tmp$batch
  condition3 <- any(is_dropped)
  expect_true(condition3)
  if(condition3){
    ##
    ## possible hemizygous deletion missing in
    ##   one component (e.g., CNP023)
    ##
    message("There are batches with fewer than 10 samples assigned to hemiz. del component with probability > 0.9")
    ##
    dropped_batches <- uniqueBatch(mod_2.3)[ is_dropped ]
    if(modelName(sb3) == "SBP3"){
      hemvar <- sigma2(sb3)[, 1]
    } else {
      hemvar <- sigma2(sb3)[, 2]
    }
    message("Find mean and variance of hemizygous deletion component")
    loc.scale.hem <- tibble(theta=theta(sb3)[2],
                            sigma2=hemvar,
                            phat=p(sb3)[2],
                            batch=dropped_batches,
                            theta.diploid=theta(mod_2.3)[dropped_batches, 2]) %>%
      mutate(delta=theta.diploid-theta)
    message("Augment data with additional hemizygous deletions")
    impdat <- impute(mod_2.3, loc.scale.hem, 1)
    obsdat <- assays(mod_2.3) %>%
      mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat,
                        impdat) %>%
      arrange(batch)
  }
  mb <- MultiBatchList(data=simdat)[[modelName(mod_2.3)]]
  mcmcParams(mb) <- mcmcParams(mod_2.3)
  theta(mb) <- theta(mod_2.3)
  theta(mb)[is_dropped, 1] <- loc.scale.hem$theta
  CNPBayes:::sigma2(mb) <- matrix(pmin(sigma2(mod_2.3)[, 1], hemvar), ncol=1)
  ##
  message("Run additional MCMC simulations on the augmented data")
  ##
  mod_2.32 <- posteriorSimulation(mb)
  mod_2.3 <- mod_2.32
  if(FALSE){
    saveRDS(mod_2.32, file="../../inst/extdata/mod_2.32.rds")
  }
  expected <- readRDS("../../inst/extdata/mod_2.32.rds")
  expect_equivalent(mod_2.32, expected)
})

test_that("augment homozygous deletions", {
  ##
  ## Do i need all 3 of these objects?
  ##
  mb.subsamp <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "mb_subsamp.rds"))
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3_1.rds"))
  mod_2.3 <- readRDS(file.path("..",
                               "..",
                               "inst",
                               "extdata",
                               "mod_2.32.rds"))
  THR <- -1
  is_pooledvar <- ncol(sigma2(mod_2.3))==1
  expect_true(is_pooledvar)
  p_ <- cbind(p(sb3)[1, 1], p(mod_2.3)) %>%
    "/"(rowSums(.))
  hdmean <- -3.887 ## computed in summarized_regions
  if(is.na(hdmean)) hdmean <- THR-1
  theta_ <- cbind(hdmean,
                  theta(mod_2.3))
  if(is_pooledvar){
    sigma2_ <- sigma2(mod_2.3)
  } else {
    sigma2_ <- cbind(sigma2(mod_2.3)[1, 2],
                     sigma2(mod_2.3))
  }
  ##
  ## Which batches have fewer than 5% subjects with
  ## homozygous deletions?
  ##
  freq.hd <- assays(mb.subsamp) %>%
    group_by(batch) %>%
    summarize(N=n(),
              n=sum(likely_hd)) %>%
    filter(n/N < 0.05)
  do_augment_homozygous <- nrow(freq.hd) > 0
  expect_true(do_augment_homozygous)
  ##
  ## Reminder that this is only evaluated when we need to augment
  ##
  if(! do_augment_homozygous ) stop()
  ##
  ## If at least one batch has fewer than 5% subjects with
  ## homozygous deletion, augment the data for homozygous deletions
  ##
  ##
  ##
  message("Augment data for homozygous deletions")
  ##
  ##
  loc.scale <- tibble(theta=hdmean,
                      sigma2=sigma2_[, 1],
                      phat=max(p(sb3)[1], 0.05),
                      batch=seq_len(nrow(theta_)))
  loc.scale <- left_join(freq.hd, loc.scale, by="batch") %>%
    select(-N)
  start.index <- length(grep("augment", id(mod_2.3))) + 1
  imp.hd <- impute(mod_2.3, loc.scale, start.index=start.index)
  condition5 <- any(isSimulated(mod_2.3))
  expect_true(condition5)
  if(condition5){
    imp.hemi <- filter(assays(mod_2.3), is_simulated)
    imp.hd <- bind_rows(imp.hd, imp.hemi)
  }
  obsdat <- assays(mb.subsamp) %>%
    mutate(is_simulated=FALSE)
  simdat <- bind_rows(obsdat, imp.hd) %>%
    arrange(batch)
  if(!exists("simdat")){
    simdat <- assays(mb.subsamp) %>%
      mutate(is_simulated=FALSE)
  }
  if(FALSE){
    saveRDS(simdat, file="../../inst/extdata/simdat_1.rds")
  }
  expected <- readRDS("../../inst/extdata/simdat_1.rds")
  expect_equivalent(simdat, expected)
})

.test_that("", {
  ##
  message("Rerun 3 batch model with augmented data")
  ##
  mbl <- MultiBatchList(data=simdat)
  model <- CNPBayes:::incrementK(mod_2.3)
  mod_1.3 <- mbl[[ model ]]
  theta(mod_1.3) <- cbind(hdmean, theta(mod_2.3))
  if(is_pooledvar){
    CNPBayes:::sigma2(mod_1.3) <- sigma2(mod_2.3)
  } else {
    CNPBayes:::sigma2(mod_1.3) <- cbind(sigma2(sb3)[1], sigma2(mod_2.3))
  }
  burnin(mod_1.3) <- 200
  iter(mod_1.3) <- 0
  mod_1.3 <- tryCatch(posteriorSimulation(mod_1.3),
                      error=function(e) NULL)
  if(is.null(mod_1.3) || is.nan(log_lik(mod_1.3))){
    mbl <- replicate(10, MultiBatchList(data=fdat)[[ model ]])
    for(j in seq_along(mbl)){
      cat(".")
      mb <- mbl[[j]]
      burnin(mb) <- 100
      iter(mb) <- 0
      mb <- posteriorSimulation(mb)
      mbl[[j]] <- mb
    }
    mod_1.3 <- mbl[[ which.max(sapply(mbl, log_lik)) ]]
  }
  if(!is.null(mod_1.3)){
    internal.count <- CNPBayes:::flags(mod_1.3)$.internal.counter
    any_dropped <- TRUE
    if(internal.count < 100){
      iter(mod_1.3) <- 1000
      mod_1.3 <- posteriorSimulation(mod_1.3)
      pz <- probz(mod_1.3) %>%
        "/"(rowSums(.)) %>%
        rowMax
      ubatch <- batch(mod_1.3)[ pz >= 0.95 & z(mod_1.3) == 2] %>%
        unique
      any_dropped <- any(!batch_labels %in% ubatch)
      ## Check that variance estimates are comparable to mod_2.3
      varratio <- sigma2(mod_1.3)/sigma2(mod_2.3)
    }  else  varratio <- Inf; internal.count <- 200; any_dropped <- TRUE
  } else {
    varratio <- Inf; internal.count <- 200; any_dropped <- TRUE
  }
  condition6 <- any(varratio > 4) || internal.count >= 100 ||
    any_dropped
  expect_true(condition6)
  if(condition6){
    ##
    ## variances much larger than expected when homozygous deletions added
    ## - allow homozygous deletions to have larger variance
    ## - fix other components
    ##
    ##
    message("MultiBatch model has higher variance than needed")
    message("Fix variances and only learn hom del component")
    ##
    model <- gsub("P", "", modelName(mod_1.3))
    mod_1.3 <- mbl[[ model ]]
    hdmean <- median(dat$oned[dat$likely_hd])
    theta(mod_1.3) <- cbind(hdmean, theta(mod_2.3))
    is_pooledvar <- ncol(sigma2(mod_1.3)) == 1
    if(!is_pooledvar){
      multibatchvar <- sigma2(mod_2.3)
      s2 <- replicate(k(mod_1.3)-1, multibatchvar, simplify=FALSE) %>%
        do.call(cbind, .)
      singlebatchvar <- sigma2(sb3)[, 1]
      foldchange_singlebatchvar <- singlebatchvar/median(multibatchvar[, 1])
      if(foldchange_singlebatchvar > 5) {
        singlebatchvar <- median(multibatchvar[, 1])*5
      }
      s2 <- cbind(singlebatchvar, s2)
    } else {
      s2 <- sigma2(mod_2.3)
    }
    CNPBayes:::sigma2(mod_1.3) <- s2
    burnin(mod_1.3) <- 200
    iter(mod_1.3) <- 1000
    mod_1.3 <- mcmc_homozygous(mod_1.3)
  }
  ##
  ## Genotype using high confidence subjects
  ##
  pz <- probz(mod_1.3) %>%
    "/"(rowSums(.))
  max_zprob <- rowMax(pz)
  is_high_conf <- max_zprob > 0.95
  mean(is_high_conf)
  gmodel <- mod_1.3[ !isSimulated(mod_1.3) &  is_high_conf ]
  keep <- !duplicated(CNPBayes:::id(gmodel))
  gmodel <- gmodel[ keep ]
  snpdat2 <- snpdat[, id(gmodel) ]
  ##identical(colnames(snpdat), id(final.model))
  clist <- CnList(gmodel)
  (stats <- baf_loglik(clist, snpdat2))
  mapping(gmodel) <- strsplit(stats$cn.model[1], ",")[[1]]

  maxpz <- probz(gmodel) %>%
    "/"(rowSums(.)) %>%
    rowMax
  bafdat <- assays(snpdat)[["baf"]] %>%
    as_tibble() %>%
    mutate(rsid=rownames(snpdat)) %>%
    gather("id", "BAF", -rsid)
  ##
  ## tibble of copy number probabilities
  ##
  cndat <- tibble(id=id(gmodel),
                  batch=batch(gmodel),
                  oned=oned(gmodel),
                  pz=maxpz,
                  z=map_z(gmodel)) %>%
    mutate(cn=mapping(gmodel)[z]) %>%
    mutate(cn=factor(cn))
  bafdat <- left_join(bafdat, cndat, by="id")
  cnlabels <- paste0(seq_len(k(gmodel)),
                     "%->%",
                     mapping(gmodel))
  labs <- as.expression(parse(text=cnlabels))
  xlab <- expression(paste("\n", "Mixture component"%->%"Copy number"))
  ##
  ##
  ##
  message("Write mapping to model file")
  mapping(mod_1.3) <- mapping(gmodel)
  expect_identical(mapping(mod_1.3), c("0", "2", "2"))
  saveRDS(mod_1.3, file=file.path(modeldir, paste0(CNP, ".rds")))
})

## Plotting code


.test_that("", {
  xlimit <- range(dat$oned)
  if(diff(xlimit) < 4){
    xlimit <- c(-3, 1)
  }
  a <- ggMixture(mb.subsamp, bins=200) +
    xlim(xlimit) +
    geom_vline(xintercept=THR, color="gray", linetype="dashed")

  xlimit <- range(oned(sb3))
  b <- ggMixture(sb3) + xlim(xlimit)
  if(FALSE) print(b)

  A <- ggMixture(mod_1.3) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n")
  ## predictive densities excluding simulated data
  A2 <- ggMixture(mod_1.3[ !isSimulated(mod_1.3) ]) +
    xlab(expression(paste("Median ", log[2], " R ratio"))) +
    ylab("Density\n")

  bafdat2 <- filter(bafdat, pz > 0.9)
  ##
  ## PLot BAFs
  ##
  B <- ggplot(bafdat2, aes(factor(z), BAF)) +
    geom_hline(yintercept=c(0, 1/3, 0.5, 2/3, 1), color="gray95") +
    geom_jitter(aes(color=pz), width=0.1, size=0.3) +
    scale_y_continuous(expand=c(0, 0.05)) +
    scale_x_discrete(breaks=seq_len(k(gmodel)),
                     labels=labs) +
    theme(panel.background=element_rect(fill="white", color="gray30"),
          legend.key=element_rect(fill="white")) +
    xlab(xlab) +
    ylab("BAF\n") +
    guides(color=guide_legend(title="Mixture\ncomponent\nprobability"))


  pdf(figname, width=14, height=8)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1, 2, widths=c(0.6, 0.4))))
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=1))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.9, "npc")))
  print(A, newpage=FALSE)
  popViewport(2)
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=2))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.6, "npc")))
  print(B, newpage=FALSE)
  dev.off()

  pdf(figname2, width=14, height=8)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1, 2, widths=c(0.6, 0.4))))
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=1))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.9, "npc")))
  print(A2, newpage=FALSE)
  popViewport(2)
  pushViewport(viewport(layout.pos.row=1,
                        layout.pos.col=2))
  pushViewport(viewport(width=unit(0.96, "npc"),
                        height=unit(0.6, "npc")))
  print(B, newpage=FALSE)
  dev.off()
  model <- mod_1.3
})
