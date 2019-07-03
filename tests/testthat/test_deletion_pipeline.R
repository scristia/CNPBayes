context("analysis of rare deletions")

## for first step, see test_summarized_region.R
.test_that <- function(nm, expr) NULL

test_that("augment data", {
  set.seed(123)
  mp <- McmcParams(iter=500, burnin=400)
  seeds <- sample(seq_len(10000), 263, replace=TRUE)
  set.seed(seeds[ 1 ])
  THR <- -1
  mb.subsamp <- readRDS("../../inst/extdata/mb_subsamp.rds")
  mean_sd <- meanSdHomDel(mb.subsamp, -1)
  rare_homozygous <- sum(oned(mb.subsamp) < THR) < 5
  expect_false(rare_homozygous)
  if(rare_homozygous){
    simdat <- augment_homozygous(mb.subsamp, mean_sd, THR)
  } else {
    simdat <- assays(mb.subsamp) %>%
      arrange(batch) %>%
      mutate(homozygousdel_mean=mean_sd[1])
  }
  if(FALSE){
    saveRDS(simdat, file="../../inst/extdata/simdat.rds")
  }
  expected <- readRDS("../../inst/extdata/simdat.rds")
  expected$homozygousdel_mean <- mean_sd[1]
  expect_equivalent(simdat, expected)
})

test_that("mcmc", {
  set.seed(5)
  simdat <- readRDS("../../inst/extdata/simdat.rds")
  message("Random starts and short warmup for SBP3 model")
  sb3 <- warmup(simdat, "SBP3", "SB3")
  set.seed(2463)
  mcmcParams(sb3) <- McmcParams(iter=400, burnin=500)
  sb3 <- posteriorSimulation(sb3)
  expect_equal(as.numeric(theta(sb3)),
               c(-3.76, -0.217, -0.0019),
               tolerance=0.01)
  finished <- stop_early(sb3)
  expect_false(finished)

  ## Since not finished, keep going
  mb <- revertToMultiBatch(sb3) 
  fdat <- filter(assays(mb), oned > -1)  %>%
    mutate(is_simulated=FALSE)
  mb <- warmup(fdat, "MBP2")
  mcmcParams(mb) <- McmcParams(iter=500, burnin=50)
  mod_2.3 <- posteriorSimulation(mb)
  ok <- ok_hemizygous(sb3)
  expect_true(ok)
  if(!ok) {
    mod_2.3 <- augment_hemizygous(mod_2.3, sb3)
    mod_2.3 <- posteriorSimulation(mod_2.3)
  }
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
  mb <- augment_hemizygous(mod_2.3, sb3)
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
  set.seed(2463)
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
  mb.subsamp <- revertToMultiBatch(sb3)
  THR <- -1
  is_pooledvar <- ncol(sigma2(mod_2.3))==1
  expect_true(is_pooledvar)
  p_ <- cbind(p(sb3)[1, 1], p(mod_2.3)) %>%
    "/"(rowSums(.))
  hdmean <- homozygousdel_mean(mb.subsamp, -1)
  if(is.na(hdmean)) hdmean <- THR-1
  ##theta_ <- cbind(hdmean, theta(mod_2.3))
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
  ## Reminder that code below is only evaluated when we need to augmen
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
                      batch=batchLevels(sb3))
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
  ##  if(!exists("simdat")){
  ##    simdat <- assays(mb.subsamp) %>%
  ##      mutate(is_simulated=FALSE)
  ##  }
  if(FALSE){
    saveRDS(simdat, file="../../inst/extdata/simdat_1.rds")
  }
  ##expected <- readRDS("../../inst/extdata/simdat_1.rds")
  ##expect_equivalent(simdat, expected)
})

test_that("Fit multi-batch model with all components", {
  set.seed(2463)
  simdat <- readRDS("../../inst/extdata/simdat_1.rds")
  mod_2.3 <- readRDS("../../inst/extdata/mod_2.3.rds")
  is_pooledvar <- ncol(sigma2(mod_2.3))==1
  hdmean <- -3.887 ## computed in summarized_regions
  batch_labels <- assays(mod_2.3)$batch_labels
  batch_labels <- batch_labels %>%
    factor(., levels=unique(.)) %>%
    as.integer(.) %>%
    unique(.) %>%
    sort()
  ##
  message("Rerun 3 batch model with augmented data")
  ##
  mbl <- MultiBatchList(data=simdat)
  model <- incrementK(mod_2.3)
  mod_1.3 <- mbl[[ model ]]
  theta(mod_1.3) <- cbind(hdmean, theta(mod_2.3))
  if(is_pooledvar){
    sigma2(mod_1.3) <- sigma2(mod_2.3)
  } else {
    sigma2(mod_1.3) <- cbind(sigma2(sb3)[1], sigma2(mod_2.3))
  }
  burnin(mod_1.3) <- 200
  iter(mod_1.3) <- 0
  mod_1.3 <- tryCatch(posteriorSimulation(mod_1.3),
                      error=function(e) NULL)
  bad_start <- is.null(mod_1.3) || is.nan(log_lik(mod_1.3))
  expect_false(bad_start)
  if(bad_start){
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
    internal.count <- flags(mod_1.3)$.internal.counter
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
    }  else  {
      varratio <- Inf; internal.count <- 200; any_dropped <- TRUE
    }
  } else {
    varratio <- Inf; internal.count <- 200; any_dropped <- TRUE
  }
  expect_equal(internal.count, 0)
  expect_false(any_dropped)
  bad_pooled_variance <- any(varratio > 4) || internal.count >= 100 ||
    any_dropped
  expect_false(bad_pooled_variance)
  if(FALSE){
    saveRDS(mod_1.3, file="../../inst/extdata/mod_1.3.rds")
  }
  expected <- readRDS("../../inst/extdata/mod_1.3.rds")
  expect_equivalent(mod_1.3, expected)
})

.test_that("only update homozygous component", {
  set.seed(2463)
  simdat <- readRDS("../../inst/extdata/simdat_1.rds")
  mod_1.3 <- readRDS("../../inst/extdata/mod_1.3.rds")
  mod_2.3 <- readRDS("../../inst/extdata/mod_2.3.rds")
  sb3 <- readRDS(file.path("..",
                           "..",
                           "inst",
                           "extdata",
                           "sb3_1.rds"))
  hdmean <- -3.887 ## computed in summarized_regions
  bad_pooled_variance <- TRUE
  mbl <- MultiBatchList(data=simdat)
  if(bad_pooled_variance){
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
    ##hdmean <- median(dat$oned[dat$likely_hd])
    theta(mod_1.3) <- cbind(hdmean, theta(mod_2.3))
    is_pooledvar <- ncol(sigma2(mod_1.3)) == 1
    expect_false(is_pooledvar)
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
  if(FALSE){
    saveRDS(mod_1.3, file="../../inst/extdata/mod_1.3_2.rds")
  }
  expected <- readRDS("../../inst/extdata/mod_1.3_2.rds")
  expect_equivalent(mod_1.3, expected)
})

test_that("genotype", {
  ##
  ## Genotype using high confidence subjects
  ##
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  i <- 1
  se <- cnp_se[i, ]
  snpdat <- snp_se[overlapsAny(snp_se, se), ]
  mod_1.3 <- readRDS("../../inst/extdata/mod_1.3.rds")
  pz <- probz(mod_1.3) %>%
    "/"(rowSums(.))
  max_zprob <- rowMax(pz)
  is_high_conf <- max_zprob > 0.95
  ##mean(is_high_conf)
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
  ##
  ##
  ##
  message("Write mapping to model file")
  mapping(mod_1.3) <- mapping(gmodel)
  expect_identical(mapping(mod_1.3), c("0", "2", "2"))
  if(FALSE){
    saveRDS(mod_1.3, file="../../inst/extdata/CNP_001.rds")
  }
  ##expected <- readRDS("~/Dropbox/Labs/klein/cnpbayes-paper/cnp.models/inst/extdata/CNP_001.rds")
  expected <- readRDS("../../inst/extdata/CNP_001.rds")
  expect_equivalent(mod_1.3, expected)
  th.expected <- theta(expected)
  ##dimnames(th.expected) <- NULL
  expect_equal(theta(mod_1.3), th.expected, tolerance=0.05)
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
  cnlabels <- paste0(seq_len(k(gmodel)),
                     "%->%",
                     mapping(gmodel))
  labs <- as.expression(parse(text=cnlabels))
  xlab <- expression(paste("\n", "Mixture component"%->%"Copy number"))
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
