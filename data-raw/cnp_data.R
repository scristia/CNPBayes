simulateFindCnpData <- function(grl, snp_exp2, cnvs){
  log.r.empirical <- lrr(snp_exp2)[, "6872.csv"]
  baf.empirical <- baf(snp_exp2)[, "6872.csv"]
  gr <- reduce(unlist(grl))
  indices.deletion <- subjectHits(findOverlaps(gr, snp_exp2))
##   Assume that all the other observations are diploid (even though
##   this is not true for this sample)
  indices.diploid <- seq_along(snp_exp2)[-indices.deletion]
  rr <- rowRanges(snp_exp2)
  nr <- length(rr)
  nc <- 35
  b.a.f <- log.r.ratios <- matrix(NA, nr, nc)
  b <- r <- rep(NA, nr)
  set.seed(123)
  for(i in seq_len(nc)){
    if(i <= 25){
      ##sample with deletion
      g.cnv <- cnvs[i]
      ##not correct!
      g.cnv <- consensusCNP(grl[1], max.width=5e6)
      cnv.region <- consensusCNP(grl[1], max.width=5e6)
      J <- subjectHits(findOverlaps(g.cnv, rr))
      i.deletion <- sample(indices.deletion, length(J), replace=TRUE)
      r[J] <- log.r.empirical[i.deletion]
      b[J] <- baf.empirical[i.deletion]
      ndiploid <- length(snp_exp2) - length(J)
      i.diploid <- sample(indices.diploid, ndiploid, replace=TRUE)
      r[-J] <- log.r.empirical[i.diploid]
      b[-J] <- baf.empirical[i.diploid]
    } else {
      ##diploid sample
      i.diploid <- sample(indices.diploid, length(r), replace=TRUE)
      r <- log.r.empirical[i.diploid]
      b <- baf.empirical[i.diploid]
    }
    b.a.f[, i] <- b
    log.r.ratios[, i] <- r
  }
  dimnames(log.r.ratios) <- dimnames(b.a.f) <- list(rownames(snp_exp2),
                                                    paste0("sample", 1:35))

  return(list(cn=log.r.ratios, baf=b.a.f))
}
library(VanillaICE)
library(CNPBayes)
data(snp_exp)
region <- GRanges("chr22", IRanges(15e6, 22e6))
hits <- findOverlaps(region, rowRanges(snp_exp))
snp_exp2 <- snp_exp[subjectHits(hits), ]
## simulate 25 samples with deletion and 10 samples without deletion
grl <- readRDS(system.file("extdata", "grl_deletions.rds", package="CNPBayes"))
lrr_baf <- simulateFindCnpData(grl, snp_exp2, grl)
sim.se <- SnpArrayExperiment(cn=lrr_baf$cn, baf=lrr_baf$baf, 
                             rowRanges=rowRanges(snp_exp2))

## depend only on SummarizedExperiment
se <- readRDS("~/Software/CNPBayes/inst/extdata/simulated_se.rds")
se2 <- SummarizedExperiment(assays=SimpleList(cn=lrr(se),
                                              baf=baf(se)),
                            rowRanges=rowRanges(se))
saveRDS(se2, file="~/Software/CNPBayes/inst/extdata/simulated_se.rds")
