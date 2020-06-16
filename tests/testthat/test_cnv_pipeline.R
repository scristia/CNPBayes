context("CNV model selection")

.test_that <- function(msg, expr) NULL

##
## Common deletion
##
test_that("common deletion", {
    path <- system.file("extdata", package="CNPBayes")
    set.seed(555)
    cnp_se <- readRDS(file.path(path, "cnp_se.rds"))["CNP_022", ]
    snp_se <- readRDS(file.path(path, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    path2 <- file.path(path, "CNP_022")
    mb.subsamp <- readRDS(file.path(path2, "mb_subsamp.rds"))
    ## batch effect is subtle.  Try both SB and MB
    ##ggMixture(mb.subsamp)
    if(FALSE){
        ggMixture(mb.subsamp, bins=300) +
            coord_cartesian(xlim=c(-1, 1)) +
            theme_bw() +
            geom_vline(xintercept=0)
        ## mode of 3rd batch is shifted to right
        mb2 <- toSingleBatch(mb.subsamp)
        ggMixture(mb2, bins=300) +
            coord_cartesian(xlim=c(-1.5, 1)) +
            theme_bw() +
            geom_vline(xintercept=0)        
    }
    ## should do more iterations in practice
    mp <- McmcParams(iter=400, burnin=100)
    model <- homdel_model(mb.subsamp, mp, Nrep=3)
    gmodel <- genotype_model(model, snp_se)
    expect_identical(mapping(gmodel), c("0", "1", "2"))
    expect_identical(modelName(gmodel), "MBP3")
})

.test_that("deletion_models", {
    library(SummarizedExperiment)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_029")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_029", ]
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    mp <- McmcParams(iter=400, burnin=500)
    set.seed(5)
    ## Only identified when I did 10 starts
    ##
    ## -- needs augmentation at duplicated allele
    ##
    model <- deletion_models(mb.subsamp, snp_semp, Nrep=10)
    if(FALSE){
        up1.4 <- readRDS("temp8.rds")
        g <- genotype_model(up1.4, snp_se)
    }
    gmodel <- genotype_model(model, snp_se)
    expect_identical(mapping(gmodel), c("0", "1", "2", "3"))
    if(FALSE){
        model2 <- select_highconfidence_samples2(model, snp_se)
        bafdat <- join_baf_oned(model2, snp_se)
        figs <- list_mixture_plots(model, bafdat)
        mixture_layout(figs, augmented=FALSE)    
    }
})

## homozygous deletion and duplication, as well as obvious batch effects
test_that("homdeldup_model", {
    library(SummarizedExperiment)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_029")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_029", ]
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    mp <- McmcParams(iter=400, burnin=500)
    set.seed(5)
    ## Only identified when I did 10 starts
    ##
    ## -- needs augmentation at duplicated allele
    ##
    model <- homdeldup_model(mb.subsamp, mp, Nrep=10)
    if(FALSE){
        up1.4 <- readRDS("temp8.rds")
        g <- genotype_model(up1.4, snp_se)
    }
    gmodel <- genotype_model(model, snp_se)
    expect_identical(mapping(gmodel), c("0", "1", "2", "3"))
    if(FALSE){
        model2 <- select_highconfidence_samples2(model, snp_se)
        bafdat <- join_baf_oned(model2, snp_se)
        figs <- list_mixture_plots(model, bafdat)
        mixture_layout(figs, augmented=FALSE)    
    }
})

.test_that("homdeldup_steps", {
    skip("Not running homdeldup_steps")
    library(SummarizedExperiment)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_029")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_029", ]
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    mp <- McmcParams(iter=400, burnin=500)
    set.seed(5)
    ## Only identified when I did 10 starts
    ## model <- homdeldup_model(mb.subsamp, mp, Nrep=10)
    mb <- mb.subsamp
    augment <- TRUE
    assays(mb) <- augment_homozygous(mb)
    sb <- warmup(assays(mb), "SBP4", "SB4", Nrep=3)
    mcmcParams(sb) <- mp
    sb <- posteriorSimulation(sb)
    expect_false(substr(modelName(sb), 3, 3) == "P")
    finished <- stop_early(sb, 0.98, 0.98)
    expect_false(finished)
    mb.subsamp <- mb
    fdat <- filter(assays(mb.subsamp), !likely_deletion)
    mb <- warmup(fdat, "MBP3", Nrep=3)
    mcmcParams(mb) <- mp
    ##message("Fitting restricted model")
    ##trace(restricted_homhemdup, browser)
    mod_2.4 <- restricted_homhemdup(mb, mb.subsamp, mp, Nrep=10)
    mod_2.4 <- readRDS("temp5.rds")
    message("Data augmentation for homozygous deletions")
    results <- list(mod_2.4=mod_2.4, sb=sb, mb.subsamp=mb.subsamp)
    saveRDS(results, file="temp6.rds")
    ##simdat2 <- augment_rarehomdel(mod_2.4, sb, mb.subsamp)
    ##simdat2 <- readRDS("temp7.rds")
    ##mod_1.4 <- hd4comp(mod_2.4, simdat2, mb.subsamp, mp)
    up1.4 <- readRDS("temp8.rds")
})

.test_that("restricted_homhemdup", {
    skip("restricted_homhemdup")
    set.seed(123)
    args <- readRDS("temp1.rds")
    mb <- args[[1]]
    mb.subsamp <- args[[2]]
    mp <- mcmcParams(mb)
    ###trace(restricted_homhemdup, browser)
    ## tmp <- restricted_homhemdup(mb, mb.subsamp, mp, Nrep=10)
    ##
    ## Fitting restricted model (no deletions)
    ##
    mod_2.4 <- posteriorSimulation(mb)
    is_flagged <- mod_2.4@flags$.internal.counter > 40
    if(!is_flagged) return(mod_2.4)
    ##
    ## Restricted MB model failed. Try SB model
    ##
    sb3 <- warmup(assays(mb), "SBP3")
    mcmcParams(sb3) <- mp
    sb3 <- posteriorSimulation(sb3)
    ## simulate from the pooled model for each batch
    full.dat <- assays(mb.subsamp)
    ##    browser()
    ##    ggMixture(mod_2.4) +
    ##        geom_vline(xintercept=theta(sb3)[, 3] - 2*sigma(sb3)[1, 1])
    thr <- theta(sb3)[, 3] - 2*sigma(sb3)[1, 1]
    args <- list(thr=thr, sb3=sb3, mod_2.4=mod_2.4, full.dat=full.dat)
    saveRDS(args, file="temp2.rds")
    ##simdat <- augment_rareduplication(sb3, mod_2.4, full.dat, thr)
    simdat <- readRDS("temp4.rds")
    mod_2.4.2 <- MultiBatchList(data=simdat)[["MBP3"]]
    simdat2 <- augment_rarehemdel(sb3,
                                  mod_2.4.2,
                                  full_data=full.dat)
    expect_true(all(tail(simdat2)$likely_deletion == FALSE))
    filtered.dat <- filter(simdat2, !likely_deletion)
    expect_identical(tail(filtered.dat), tail(simdat2))
    mb <- warmup(filtered.dat, "MBP3")
    mcmcParams(mb) <- mp
    mod_2.4 <- posteriorSimulation(mb)
    saveRDS(mod_2.4, "temp5.rds")
})

.test_that("augment_rarehomdel", {
    skip("augment_rarehomdel")
    args <- readRDS("temp6.rds")
    restricted <- args[["mod_2.4"]]
    sb4 <- args[["sb"]]
    mb.subsamp <- args[["mb.subsamp"]]

    p_ <- cbind(p(sb4)[1, 1], p(restricted)) %>%
        "/"(rowSums(.))
    dat <- assays(mb.subsamp)
    hdmean <- median(dat$oned[dat$likely_deletion])
    hdvar <- var(dat$oned[dat$likely_deletion])
    if(is.na(hdmean)) hdmean <- -4
    theta_ <- cbind(hdmean, theta(restricted))
    is_pooledvar <- ncol(sigma(restricted)) == 1
    expect_true(is_pooledvar)
    sigma2_ <- sigma2(restricted)
    freq.hd <- assays(mb.subsamp) %>%
        group_by(batch) %>%
        summarize(N=n(),
                  n=sum(likely_deletion)) %>%
        filter(n/N < 0.05)
    expect_true(nrow(freq.hd) > 0)
    loc.scale <- tibble(theta=hdmean,
                        sigma2=sigma2_[, 1],
                        phat=max(p(sb4)[1], 0.05),
                        batch=seq_len(nrow(theta_)))
    loc.scale <- left_join(freq.hd, loc.scale, by="batch") %>%
        select(-N)
    start.index <- length(grep("augment", id(restricted))) + 1
    imp.hd <- impute(restricted, loc.scale, start.index=start.index)
    expect_true(any(isSimulated(restricted)))
    imp1 <- filter(assays(restricted), is_simulated)
    imp.hd <- bind_rows(imp.hd, imp1)
    obsdat <- assays(mb.subsamp) %>%
        mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat, imp.hd) %>%
        arrange(batch)
    saveRDS(simdat, "temp7.rds")
})

.test_that("hd4comp", {
    skip("hd4comp")
    simdat2 <- readRDS("temp7.rds")
    args <- readRDS("temp6.rds")
    mod_2.4 <- args[["mod_2.4"]]
    mb.subsamp <- args[["mb.subsamp"]]
    mp <- mcmcParams(mod_2.4)
    
    model <- incrementK(mod_2.4) %>%
        gsub("P", "", .)
    ##THR <- summaries(mb.subsamp)$deletion_cutoff
    THR <- deletion_midpoint(mb.subsamp)
    mod_1.4 <- MultiBatchList(data=simdat2)[[ model ]]
    hdmean <- homozygousdel_mean(mb.subsamp, THR)
    hdvar <- homozygousdel_var(mb.subsamp, THR)
    theta(mod_1.4) <- cbind(hdmean, theta(mod_2.4))
    V <- matrix(sigma2(mod_2.4)[, 1],
                nrow(sigma2(mod_2.4)), 3,
                byrow=FALSE)
    sigma2(mod_1.4) <- cbind(hdvar, V)
    mcmcParams(mod_1.4) <- mp
    ##mod_1.4 <- mcmc_homozygous(mod_1.4)
    up1.4 <- posteriorSimulation(mod_1.4)
    up1.4
    saveRDS(up1.4, "temp8.rds")
})

.test_that("compute_density", {
    skip("compute_density")
    args <- readRDS("temp2.rds")
    thr <- args[["thr"]]
    mb <- args[["mod_2.4"]]
    batches <- batchLevels(mb)
    densities <- vector("list", length(batches))
    for(i in seq_along(batches)){
        ##
        ## Coarse
        ##
        index <- which(batch(mb)==i & oned(mb) > thr)
        if(length(index) > 2){
            tmp <- density(oned(mb)[index]) %>%
                "["(1:2)
            tmp <- as_tibble(tmp)
        } else tmp <- NULL
        ##
        ## Fine
        ##
        tmp2 <- density(oned(mb)[batch(mb)==i & oned(mb) < thr]) %>%
            "["(1:2) %>%
            as_tibble()
        densities[[i]] <- list(coarse=as_tibble(tmp),
                               fine=as_tibble(tmp2))
    }
    densities    
    saveRDS(densities, "temp3.rds")
})

.test_that("augment_rareduplication", {
    skip("augment_rareduplication")
    args <- readRDS("temp2.rds")
    sb3 <- args[["sb3"]]
    mod_2.4 <- args[["mod_2.4"]]
    full_data <- args[["full.dat"]]
    THR <- args[["thr"]]
    
    densities <- readRDS("temp3.rds")
    modes <- round(compute_modes(densities), 3)
    tab <- tibble(modes_=modes, batch=paste0("Batch ", seq_along(modes)))
    fig <- ggMixture(mod_2.4) +
        geom_vline(data=tab, aes(xintercept=modes_))
    
    loc.scale <- tibble(theta=theta(sb3)[3],
                        sigma2=sigma2(sb3),
                        phat=max(p(sb3)[1, 3], 0.05),
                        batch=seq_len(numBatch(mod_2.4)))
    loc.scale$theta <- loc.scale$theta + modes
    start.index <- length(grep("augment", id(mod_2.4))) + 1
    imp.dup <- impute(mod_2.4, loc.scale,
                      start.index=start.index) %>%
        mutate(likely_deletion=FALSE)
    obsdat <- full_data %>%
        mutate(is_simulated=FALSE)
    simdat <- bind_rows(obsdat, imp.dup) %>%
        arrange(batch)
    mb <- MultiBatch(data=simdat)
    locs <- loc.scale %>%
        mutate(batch=paste0("Batch ", batch))
    fig <- ggMixture(mb) +
        geom_vline(data=locs, aes(xintercept=theta))
    saveRDS(simdat, file="temp4.rds")
})

##
## hemizygous deletion. possibly duplication
##
test_that("hemideletion_model", {
    set.seed(9209)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_014")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se["CNP_014"])
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    ##  hemizgyous deletion and possibly duplication
    ##  slight batch effect
    ##  ggMixture(mb.subsamp) + xlim(c(-3, 1)) + geom_vline(xintercept=-0.25)
    ## summaries(mb.subsamp)$deletion_cutoff <- -0.25
    assays(mb.subsamp) <- assays(mb.subsamp) %>%
        mutate(likely_deletion = oned <= -0.25)
    mp <- McmcParams(iter=400, burnin=500)
    mb1 <- hemdel_model(mb.subsamp, mp, Nrep=3)
    expect_identical(modelName(mb1), "MBP2")
    ## The above fits the data ok, but we might want to explore
    ## duplication models as well. This CNP provides motivation for
    ## doing so.
    ##
    ## Check for duplication
    ##
    ## we mess up the first batch -- variance is too big
    mb2 <- hemdeldup_model2(mb.subsamp, mp, Nrep=10)
    expect_identical(k(mb2), 3L)
    
    model.list <- list(mb1, mb2)
    mb1 <- genotype_model(mb1, snp_se)
    mb2 <- genotype_model(mb2, snp_se)
    model.list <- list(mb1, mb2)
    model <- choose_model(model.list)
    expect_identical(modelName(model), "MBP2")
    ## ggMixture(model) + xlim(c(-3, 1))
    ## tmp <- model.list[[2]][!isSimulated(model.list[[2]])]
    ## ggMixture(mb.subsamp) + xlim(c(-3, 1))
    ## ggMixture(tmp) + xlim(c(-3, 1))
    freq <- as.integer(table(map_z(model)))
    expect_true(gap::hwe(c(0, freq), data.type="count")$p.x2 < 0.05)
})

test_that("Hemizygous deletion +/- duplication pipeline", {
    ##
    ## For hemizygous deletions/duplications, cutoff should depend on batch -- assume it is the same distance from the mode
    ##   - easy to calculate mode of diploid component across batches
    ##   - from single batch model, we have an overall duplication mode that could be shifted by batch
    ##
    set.seed(9209)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_147")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_147", ]
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    ## obvious batch effect. Possibly 1 hemizygous deletion
    ## ggMixture(mb.subsamp) + geom_vline(xintercept=-0.25)
    ##summaries(mb.subsamp)$deletion_cutoff <- -0.25
    assays(mb.subsamp) <- assays(mb.subsamp) %>%
        mutate(likely_deletion=oned < -0.25)
    mp <- McmcParams(iter=400, burnin=500)
    set.seed(9209)
    mb1 <- hemdel_model(mb.subsamp, mp, Nrep=3)
    mb2 <- hemdeldup_model2(mb.subsamp, mp, Nrep=3)
    mb1 <- genotype_model(mb1, snp_se)
    mb2 <- genotype_model(mb2, snp_se)
    model.list <- list(mb1, mb2)
    model <- choose_model(model.list)
    expect_identical(mapping(model), c("2", "3"))
    expect_identical(modelName(model), "MBP2")
})


