context("BAF likelihood")

.test_that <- function(expr, name) NULL

.test_that("baf likelihood", {
  library(tidyverse)
  library(purrr)
  library(magrittr)
  library(SummarizedExperiment)
  data(SingleBatchModelExample, package="CNPBayes")
  sb <- SingleBatchModelExample
  if(FALSE){
    set.seed(123)
    sb <- posteriorSimulation(sb)
    SingleBatchModelExample <- sb
    save(SingleBatchModelExample, file="CNPBayes/data/SingleBatchModelExample.rda", compression_level=9)
  }
  ##
  ## Assume components are 0, 1, 2 and that the CNP region contains 2 SNPs
  true.model <- CopyNumberModel(sb)
  mapping(true.model) <- c("0", "1", "2")
  ##
  ## Simulate bafs for these samples
  ##
  N <- length(y(sb))
  simulateGt <- function(n, p.aa, p.ab, p.bb, gt.model){
    cn <- copyNumber(gt.model)
    p <- c(p.aa, p.ab, p.bb)
    freq <- table(cn)
    gt <- rep(NA, length(cn))
    gt [ cn == 2 ] <- sample(1:3, freq[3], prob=p, replace=TRUE)
    gt [ cn == 1 ] <- sample(c(1, 3), freq[2], prob=p[c(1, 3)],
                             replace=TRUE)
    gt [ cn == 0 ] <- sample(1:3, freq[1], prob=c(1, 1, 1),
                             replace=TRUE)
    list(gt)
  }
  simulateBaf <- function(gt, gt.model){
    cn <- copyNumber(gt.model)
    params <- shapeParams()
    gt <- gt[[1]]
    b <- rep(NA, length(gt))
    ## AA
    freq <- table(gt)
    aa <- params["zeroB", ]
    ab <- params["balanced", ]
    bb <- params["allB", ]
    n0 <- sum(cn == "0")
    browser()
    b[ gt == 1 ] <- rbeta(freq[1], aa[1], aa[2])
    if(any(gt == 2))
      b[ gt == 2 ] <- rbeta(freq[2], ab[1], ab[2])
    if(any(gt == 3))
      b[ gt == 3 ] <- rbeta(freq[3], bb[1], bb[2])
    ## overwrite BAFs for homozygous deletion
    b[ cn == "0" ] <- rbeta(n0, 1, 1)
    list(b)
  }
  ## simulate data under true model
  tab <- tibble(snp=rep(1:4, N),
                p.b=rep(c(0.1, 0.4, 0.8, 0.95), N),
                p.a=1-p.b,
                p.bb=p.b^2,
                p.aa=p.a^2,
                p.ab=2*p.a*p.b) %>%
    group_by(snp) %>%
    summarize(n=n(),
              p.b=unique(p.b),
              p.aa=unique(p.aa),
              p.ab=unique(p.ab),
              p.bb=unique(p.bb),
              gt=simulateGt(n, p.aa, p.ab, p.bb, true.model),
              baf=simulateBaf(gt, true.model))
  id <- paste0("sample", seq_along(y(sb)))
  ## B is returned by modelProb
  BAF.matrix <- do.call(cbind, tab$baf) %>% t %>%
    as.matrix %>%
    set_colnames(id)
  GT.matrix <- do.call(cbind, tab$gt) %>% t %>%
    as.matrix %>%
    set_colnames(id)
  ## make up fake granges
  snpdat <- SummarizedExperiment(assays=SimpleList(baf=BAF.matrix,
                                                   GT=GT.matrix))
  gt <- genotypes(snpdat)
  ## the empirical probabilities will not be the same as the true values since this is calculated over all samples, including samples with homozygous and hemizygous deletions
  p.b <- p_b(gt)
  ## calculate the likelihood for the true model
  model1 <- CopyNumberModel(sb)
  mapping(model1) <- c("0", "1", "2")
  pwr <- 1/10
  pz <- probCopyNumber(model1) %>%
    as.tibble %>%
    set_colnames(paste0("prob_cn", unique(mapping(model1)))) %>%
    mutate(id=id)
  cols <- paste0("prob_cn", unique(mapping(model1)))
  B <- do.call(cbind, tab$baf) %>% t %>%
    as.tibble %>%
    set_colnames(id) %>%
    mutate(p.b=tab$p.b) %>%
    gather("id", "baf", -p.b) %>%
    mutate(p0=d0(baf),
           p1=d1(baf, p.b),
           p2=d2(baf, p.b),
           p3=d3(baf, p.b),
           p4=d4(baf, p.b)) %>% ## end mixture probs
    group_by(id) %>%
    summarize(baf=mean(baf, na.rm=TRUE),
              p0=prod(p0, na.rm=TRUE)^pwr,
              p1=prod(p1, na.rm=TRUE)^pwr,
              p2=prod(p2, na.rm=TRUE)^pwr,
              p3=prod(p3, na.rm=TRUE)^pwr,
              p4=prod(p4, na.rm=TRUE)^pwr) %>%
    left_join(pz, by="id") %>%
    ## from .pBAF_012 
    select(c("baf", "p0", "p1", "p2", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn1*p1 + prob_cn2*p2)
  ll.true <- sum(log(B$prob))
  ##
  ## Calculate log likelihood for model where p.b is estimated empirically
  ##
  B2 <- do.call(cbind, tab$baf) %>% t %>%
    as.tibble %>%
    set_colnames(id) %>%
    mutate(p.b=p_b(gt),
           snp.id=1:4) %>%
    gather("id", "baf", -c(p.b, snp.id)) %>%
    mutate(p0=d0(baf),
           p1=d1(baf, p.b),
           p2=d2(baf, p.b),
           p3=d3(baf, p.b),
           p4=d4(baf, p.b)) %>% ## end mixture probs
    group_by(id) %>%
    summarize(baf=mean(baf, na.rm=TRUE),
              p0=prod(p0, na.rm=TRUE)^pwr, ## p0 is the product across SNPs for one sample
              p1=prod(p1, na.rm=TRUE)^pwr,
              p2=prod(p2, na.rm=TRUE)^pwr,
              p3=prod(p3, na.rm=TRUE)^pwr,
              p4=prod(p4, na.rm=TRUE)^pwr) %>%
    left_join(pz, by="id") %>%
    ## from .pBAF_012 
    select(c("baf", "p0", "p1", "p2", cols)) %>%
    mutate(prob=prob_cn0*p0 + prob_cn1*p1 + prob_cn2*p2)
  expect_false(identical(B$p1, B2$p1))
  ll.model1 <- sum(log(B2$prob))
  ## because p.b can be a poor estimate of the true allele frequency
  ## when deletion is common
  expect_true(ll.model1 < ll.true)
  tmp <- .modelProb(model1, snpdat)
  expect_identical(tmp$baf, B2$baf)
  expect_identical(tmp$p1, B2$p1)
  expect_identical(tmp$p2, B2$p2)
  ll <- modelProb(model1, snpdat)
  expect_identical(ll, ll.model1)
})
