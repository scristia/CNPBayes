context("CNV model selection")

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
    ##ggMixture(mb.subsamp)
    ## should do more iterations in practice
    mp <- McmcParams(iter=400, burnin=100)
    ## model.list <- deletion_models(mb.subsamp, snp_se, mp)
    ## model <- choose_model(model.list, mb.subsamp)
    model <- homdel_model(mb.subsamp, mp)
    expect_identical(mapping(model), c("0", "1", "2"))
    expect_identical(modelName(model), "MBP3")
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
    ## obvious batch effect:
    ## ggMixture(mb.subsamp)
    ## since obvious batch effect, we could skip SB models
    model <- homdeldup_model(mb.subsamp, mp, skip_SB=TRUE)
    gmodel <- genotype_model(model, snp_se)
    expect_identical(mapping(gmodel), c("0", "1", "2", "3"))
    if(FALSE){
      model2 <- select_highconfidence_samples(model, snp_se)
      bafdat <- join_baf_oned(model, snp_se)
      figs <- list_mixture_plots(model, bafdat)
      mixture_layout(figs, augmented=TRUE)    
    }
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
    summaries(mb.subsamp)$deletion_cutoff <- -0.25
    mp <- McmcParams(iter=400, burnin=500)
    model.list <- hemideletion_models(mb.subsamp, snp_se, mp)
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
    summaries(mb.subsamp)$deletion_cutoff <- -0.25
    mp <- McmcParams(iter=400, burnin=500)
    set.seed(9209)
    model.list <- hemideletion_models(mb.subsamp, snp_se, mp, skip_SB=TRUE)
    model <- choose_model(model.list)
    expect_identical(mapping(model), c("2", "3"))
    expect_identical(modelName(model), "MBP2")
})
