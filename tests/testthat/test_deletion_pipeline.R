context("Deletion pipeline")

## CNP_001 is a somewhat rare homozygous deletion.  Hemizygous deletions overlap the diploid component, and so we undercount the number of hemizygous deletions 

prepareData <- function(){
    ##data(cnp_se, package="panc.data")
    ##data(snp_se, package="panc.data")
    extdir <- system.file("extdata", package="CNPBayes")
    cnp_se <- readRDS(file.path(extdir, "cnp_se.rds"))
    snp_se <- readRDS(file.path(extdir, "snp_se.rds"))
    g <- GRanges("chr1", IRanges(1627805, 1673809),
                 seqinfo=Seqinfo("chr1", seqlengths=249250621, genome="hg19"))
    snp_se <- snp_se[overlapsAny(snp_se, g), ]
    i <- 1
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
                  oned=assays(se)[[1]][1, ],
                  provisional_batch=colData(se)$Sample.Plate) %>%
        mutate(likely_deletion = oned < THR)
}

addBatchLabels <- function(dat, mb.subsamp){
    batches <- assays(mb.subsamp) %>%
        select(provisional_batch, batch) %>%
        group_by(provisional_batch) %>%
        summarize(batch=unique(batch), .groups='drop')
    pr.batch <- assays(mb.subsamp)$provisional_batch
    stopifnot(all(pr.batch %in% dat$provisional_batch))
    message("Adding the estimated batches to the original data")
    dat <- left_join(dat, batches, by="provisional_batch")
    dat
}

numberLowProb <- function(sb){
    pz <- probz(sb)
    maxprob <- rowMax(pz)
    number_lowprob <- sum(maxprob < 0.95)
    number_lowprob
}

test_that("homdel_model", {
    library(SummarizedExperiment)
    library(ggplot2)
    library(grid)
    library(panc.data)
    library(dplyr)
    library(readxl)
    set.seed(123)
    seeds <- sample(seq_len(10000), nrow(cnp_se), replace=TRUE)
    set.seed(seeds[ 1 ])    
    THR <- -1
    dat <- prepareData()
    xlimit <- range(dat$oned)
    if(diff(xlimit) < 4){
        xlimit <- c(-3, 1)
    }
    dat.nohd <- filter(dat, !likely_deletion)
    ##
    ## Group chemistry plates, excluding homozygous deletions
    ##
    ix <- sample(seq_len(nrow(dat.nohd)), 1000, replace=TRUE)
    message("Downsampling non-homozygous deletions")
    if(FALSE){
        mb.subsamp <- dat.nohd[ix, ] %>%
            bind_rows(filter(dat, likely_hd)) %>%
            mutate(is_simulated=FALSE) %>%
            MultiBatch("MB3", data=.) %>%
            findSurrogates(0.001, THR)
    }
    extdir <- system.file("extdata", package="CNPBayes")
    path1 <- file.path(extdir, "CNP_001")
    mb.subsamp <- readRDS(file.path(path1, "mb_subsamp.rds"))            
    expect_identical(numBatch(mb.subsamp), 5L)
    dat2 <- addBatchLabels(dat, mb.subsamp)
    if(FALSE){
        a <- ggMixture(mb.subsamp, bins=200) +
            xlim(xlimit) +
            geom_vline(xintercept=THR, color="gray", linetype="dashed")
    }
    ##
    message("Check batch size")
    ##
    batchfreq <- assays(mb.subsamp) %>%
        group_by(batch) %>%
        summarize(n=n(), .groups='drop')
    expect_false(any(batchfreq$n < 50))
    ##
    ## With a bad choice of start, we can get stuck in a local mode
    ## Look at several starts and select the most promising
    ##
    ## - 10 random starts for both SB3 and SBP3 models
    ## - additional MCMC simulations for model with highest log lik
    ##
    expect_false(sum(oned(mb.subsamp) < THR) < 5)
    simdat <- assays(mb.subsamp)
    ##
    sb3 <- warmup(simdat, "SBP3", "SB3",
                  model2.penalty=50)
    expect_identical(modelName(sb3), "SB3")
    ##
    message("Run 500 burnin and 400 additional simulations.")
    ##
    iter(sb3) <- 400
    burnin(sb3) <- 500
    sb3 <- posteriorSimulation(sb3)
    nl <- numberLowProb(sb3)
    expect_true(nl > 40)
    xlimit <- range(oned(sb3))
    b <- ggMixture(sb3) + xlim(xlimit)
    if(FALSE) print(b)
    finished <- stop_early(sb3)
    expect_false(finished)
    if(FALSE){
        rdsfile <- file.path("..", "..", "inst", "extdata",
                             "sb3_001.rds")
        saveRDS(sb3, file=rdsfile)
    }
})

test_that("fit_restricted", {
    ## Focus on non-homozygous component
    set.seed(13)
    extdir <- system.file("extdata", package="CNPBayes")
    THR <- -1
    sb <- readRDS(file.path(extdir, "sb3_001.rds"))
    mb <- revertToMultiBatch(sb)
    restricted <- fit_restricted(mb, sb, THR, model="MBP2")
    if(FALSE){
        rdsfile <- file.path("..", "..", "inst", "extdata",
                             "mb2_001.rds")
        saveRDS(restricted, file=rdsfile)
    }
    expect_equal(sum(assays(restricted)$is_simulated), 78)
    pz <- probz(restricted)
    calls <- ifelse(pz[, 1] > 0.5, 1, 2)
    nhom <- sum(oned(sb) < THR)
    freq <- c(nhom, as.integer(table(calls)))
    p <- gap::hwe(freq, data.type="count")$p.x2
    ## consistent with HWE
    expect_true(p > 0.05)
})

test_that("explore_multibatch", {
    set.seed(15)
    extdir <- system.file("extdata", package="CNPBayes")    
    sb <- readRDS(file.path(extdir, "sb3_001.rds"))
    mb <- revertToMultiBatch(sb)    
    restricted.mb <- readRDS(file.path(extdir, "mb2_001.rds"))

    hdmean <- homozygousdel_mean(sb, -1)
    mn_sd <- list(hdmean, sdRestrictedModel(restricted.mb))
    mb.observed <- mb[ !isSimulated(mb) ]
    rename <- dplyr::rename
    assays(mb.observed) <- assays(mb.observed) %>%
        rename(likely_deletion=likely_hd)
        
    model <- "MB3"
    mb1 <- filter(assays(restricted.mb),
                  isSimulated(restricted.mb)) %>%
        bind_rows(assays(mb.observed)) %>%
        arrange(batch) %>%
        MultiBatchList(data=.) %>%
        "[["(model)
    ##
    ## need to define likely_deletion in mb1 object
    ##
    simdat <- .augment_homozygous(mb1, mn_sd, -1,
                                  phat=max(p(sb)[1], 0.05))  
    full <- .mcmcWithHomDel(simdat, restricted.mb)
    f <- ggMixture(full)
    ok <- ok_model(full, restricted.mb)
    expect_true(ok)
    if(FALSE){
        saveRDS(full, file=file.path("..", "..",
                                     "inst", "extdata",
                                     "full_001.rds"))
    }
        
})

test_that("genotype_model", {
    set.seed(20)
    extdir <- system.file("extdata", package="CNPBayes")        
    snp_se <- readRDS(file.path(extdir, "snp_se.rds"))
    cnp_se <- readRDS(file.path(extdir, "cnp_se.rds"))    
    snp_se <- subsetByOverlaps(snp_se, cnp_se[1, ])
    full <- readRDS(file.path(extdir, "full_001.rds"))
    
    model <- full
    snpdat <- snp_se
    ## this seems pretty stringent
    gmodel <- select_highconfidence_samples(model, snpdat)
    expect_true(all(mapping(gmodel)=="?"))
    gmodel2 <- select_highconfidence_samples2(model, snpdat)    
    gmodel <- gmodel2
    keep <- !duplicated(id(gmodel))
    gmodel <- gmodel[ keep ]
    snpdat2 <- snpdat[, id(gmodel) ]
    ##identical(colnames(snpdat), id(final.model))
    clist <- CnList(gmodel)
    stats <- baf_loglik(clist, snpdat2)
    mapping(gmodel) <- strsplit(stats$cn.model[1], ",")[[1]]
    mapping(model) <- mapping(gmodel)
    summaries(model)$baf_loglik <- stats
    if(FALSE){
        bafdat <- join_baf_oned(model, snpdat2)
        figs <- list_mixture_plots(model, bafdat)
        mixture_layout(figs, augmented=FALSE)    
    }
    model2 <- model[!isSimulated(model)]
    freq <- as.integer(table(map_z(model2)))
    ## we are undercounting the number of hemizygous deletions
    expect_true(gap::hwe(freq, data.type="count")$p.x2 < 0.05)


    extdir <- system.file("extdata", package="CNPBayes")
    path1 <- file.path(extdir, "CNP_001")
    mb.subsamp <- readRDS(file.path(path1, "mb_subsamp.rds"))
    cnp_se <- readRDS(file.path(extdir, "cnp_se.rds"))
    snp_se <- readRDS(file.path(extdir, "snp_se.rds"))
    g <- GRanges("chr1", IRanges(1627805, 1673809),
                 seqinfo=Seqinfo("chr1", seqlengths=249250621, genome="hg19"))
    snp_se <- snp_se[overlapsAny(snp_se, g), ]    
    model3 <- cnv_models(mb.subsamp, rowRanges(cnp_se)[1], snp_se)
    expect_identical(mapping(model3), mapping(model2))
    if(FALSE){
        ## verify probabilities well calibrated with constrained mcmc
        model <- model3
        maxpz <- tibble(prob=rowMaxs(probz(model)),
                        batch=batch(model),
                        z=map_z(model))
        cutoffs <- group_by(maxpz, z) %>%
            summarize(upper_quartile=quantile(prob, 0.75))
        ggplot(maxpz, aes(z, prob)) +
            geom_jitter(width=0.1) +
            facet_wrap(~batch)
        bafdat <- join_baf_oned(model, snp_se)
        figs <- list_mixture_plots(model, bafdat)
        mixture_layout(figs, augmented=TRUE)
        mixture_layout(figs, augmented=FALSE)
    }    
})



