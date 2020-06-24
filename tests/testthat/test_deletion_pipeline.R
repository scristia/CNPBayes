context("Deletion pipeline")

## CNP_001 is a somewhat rare homozygous deletion.  Hemizygous deletions overlap the diploid component, and so we undercount the number of hemizygous deletions 

numberLowProb <- function(sb){
    pz <- probz(sb)
    maxprob <- rowMax(pz)
    number_lowprob <- sum(maxprob < 0.95)
    number_lowprob
}

test_that("sb_model", {
    library(SummarizedExperiment)
    library(ggplot2)
    library(grid)
    library(dplyr)
    library(readxl)
    set.seed(123)
    extdir <- system.file("extdata", package="CNPBayes")
    fname <- file.path(extdir, "cnp_se.rds")
    cnp_se <- readRDS(fname)
    seeds <- sample(seq_len(10000), nrow(cnp_se), replace=TRUE)
    set.seed(seeds[ 1 ])    
    fname <- file.path(extdir, "CNP_001",
                       "batched_data.rds")
    batched.data <- readRDS(fname)
    message("Check batch size")
    batchfreq <- batched.data %>%
        group_by(batch) %>%
        tally()
    expect_false(any(batchfreq$n < 50))
    ##
    ## With a bad choice of start, we can get stuck in a local mode
    ## Look at several starts and select the most promising
    ##
    ## - 5 random starts for both SB3 and SBP3 models
    ## - additional MCMC simulations for model with highest log lik
    ##
    expect_identical(sum(batched.data$likely_deletion), 10L)
    sb3 <- warmup(batched.data,
                  "SBP3", "SB3",
                  model2.penalty=50,
                  Nrep=5)
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

test_that("explore_multibatch", {
    ## Focus on non-homozygous component
    set.seed(13)
    extdir <- system.file("extdata", package="CNPBayes")
    sb <- readRDS(file.path(extdir, "sb3_001.rds"))
    mb <- revertToMultiBatch(sb)
    expect_identical(iter(mb), iter(sb))
    expect_identical(burnin(mb), burnin(sb))        
    mbr <- assays(mb) %>%
        filter(!likely_deletion) %>%
        MultiBatch(data=.)
    mcmcParams(mbr) <- mcmcParams(mb)
    ok <- ok_hemizygous(sb)
    expect_false(ok)
    mbr <- augment_rarecomponent(mbr, sb)
    restricted <- fit_restricted2(mbr, model="MBP2")
    if(FALSE){
        ggMixture(restricted)
        rdsfile <- file.path("..", "..", "inst", "extdata",
                             "mb2_001.rds")
        saveRDS(restricted, file=rdsfile)
    }
})

test_that("homdel_model", {
    extdir <- system.file("extdata", package="CNPBayes")
    fname <- file.path(extdir, "CNP_001",
                       "batched_data.rds")
    batched.data <- readRDS(fname)
    mb <- MultiBatch(data=batched.data)
    ## Start of homdel_model
    assays(mb) <- augment_homozygous(mb)
    mp <- McmcParams(iter=500, burnin=100)
    sb3 <- evaluate_sb3(mb, mp)
    expect_false(stop_early(sb3) || numBatch(mb) == 1)
    final <- explore_multibatch(sb3)
    f <- ggMixture(final)
    expect_identical(modelName(final), "MBP3")
    if(FALSE){
        saveRDS(final, file=file.path("..", "..",
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
    gmodel <- select_highconfidence_samples2(model, snpdat)    
    expect_false(any(duplicated(id(gmodel))))
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
    if(FALSE){
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
    }
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



