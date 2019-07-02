context("deletion pipeline")

test_that("summarize region", {
  library(SummarizedExperiment)
  library(panc.data)
  library(CNPBayes)
  library(readxl)
  ##
  ## summarize_region code chunk
  ##
  data(cnp_se, package="panc.data")
  data(snp_se, package="panc.data")
  g <- GRanges("chr1", IRanges(1627805, 1673809),
               seqinfo=Seqinfo("chr1", 249250621,
                               genome="hg19"))
  snp_se <- snp_se[overlapsAny(snp_se, g), ]
  i <- 1
  set.seed(123)
  seeds <- sample(seq_len(10000), nrow(cnp_se), replace=TRUE)
  set.seed(seeds[ i ])
  se <- cnp_se[i, ]
  CNP <- rownames(se)
  basedir <- tempdir()
  figname <- file.path(basedir, paste0(CNP, ".pdf"))
  modeldir <- file.path(basedir, "cnp.models/inst/extdata")
  if(!dir.exists(modeldir)) dir.create(modeldir, recursive=TRUE)
  colnames(cnp_se) <- colnames(snp_se)
  snpdat <- snp_se[overlapsAny(snp_se, se), ]
  ##
  ## Flag homozygous deletions
  ##
  message("Flagging apparent homozygous deletions")
  THR <- -1
  dat <- tibble(id=colnames(se),
                oned=assays(se)[["MEDIAN"]][1, ],
                provisional_batch=colData(se)$Sample.Plate) %>%
    mutate(likely_hd = oned < THR)
  dat.nohd <- filter(dat, !likely_hd)
  ##
  ## Group chemistry plates, excluding homozygous deletions
  ##
  ix <- sample(seq_len(nrow(dat.nohd)), 1000, replace=TRUE)
  message("Downsampling non-homozygous deletions")
  mb.subsamp <- dat.nohd[ix, ] %>%
    bind_rows(filter(dat, likely_hd)) %>%
    mutate(is_simulated=FALSE) %>%
    MultiBatch("MB3", data=.) %>%
    findSurrogates(0.001, THR)

  ##print(a)
  batches <- assays(mb.subsamp) %>%
    group_by(provisional_batch) %>%
    summarize(batch=unique(batch))
  pr.batch <- assays(mb.subsamp)$provisional_batch
  stopifnot(all(pr.batch %in% dat$provisional_batch))
  dat <- left_join(dat, batches, by="provisional_batch")
  ##
  ## We need the number in `hdmean` later. Where to keep it?
  ##
  hdmean <- median(dat$oned[dat$likely_hd])
  expect_equal(hdmean, -3.887, tolerance=0.001)
  ##
  message("Check batches")
  ##
  batchfreq <- assays(mb.subsamp) %>%
    group_by(batch) %>%
    summarize(n=n())
  if(any(batchfreq$n < 50)){
    batchfreq <- filter(batchfreq, n < 50)
    adat <- assays(mb.subsamp) %>%
      filter(!batch %in% batchfreq$batch)
    bdat <- filter(dat, batch %in% batchfreq$batch)
    adat2 <- bind_rows(adat, bdat)
    mb.subsamp <- MultiBatch("MB2", data=adat2)
  }

  if(FALSE){
    saveRDS(mb.subsamp, file="../../inst/extdata/mb_subsamp.rds")
    ##saveRDS(dat, file="../../inst/extdata/dat.rds")
  }
  expected <- readRDS("../../inst/extdata/mb_subsamp.rds")
  expect_equivalent(mb.subsamp, expected)
})
