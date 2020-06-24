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
    expect_true(!any(isSimulated(gmodel)))

    ## this calls 'deletion_models'
    model2 <- cnv_models(mb.subsamp,
                         rowRanges(cnp_se)["CNP_022"],
                         snp_se)
    expect_identical(nrow(model2), nrow(gmodel))
    z_gmodel <- map_z(gmodel)
    z_model2 <- map_z(model2)
    expect_identical(modelName(gmodel), modelName(model2))
    expect_true(mean(z_gmodel == z_model2) > 0.999)
})

test_that("deletion_models", {
    library(SummarizedExperiment)
    path1 <- system.file("extdata", package="CNPBayes")
    path <- file.path(path1, "CNP_029")
    cnp_se <- readRDS(file.path(path1, "cnp_se.rds"))["CNP_029", ]
    snp_se <- readRDS(file.path(path1, "snp_se.rds"))
    snp_se <- subsetByOverlaps(snp_se, cnp_se)
    mb.subsamp <- readRDS(file.path(path, "mb_subsamp.rds"))
    mp <- McmcParams(iter=1000, burnin=200, nStarts=4)
    set.seed(5)
    model2 <- homdeldup_model(mb.subsamp, mp=mp, Nrep=10)
    ## Only identified when I did 10 starts
    ##
    ## -- needs augmentation at duplicated allele
    ##
    model.list <- deletion_models(mb.subsamp, snp_se, mp=mp)
    gmodel <- choose_model(model.list, mb.subsamp)
    expect_identical(mapping(gmodel), c("0", "1", "2", "3"))
    if(FALSE){
        model2 <- CNPBayes:::select_highconfidence_samples2(gmodel,
                                                            snp_se)
        bafdat <- join_baf_oned(gmodel, snp_se)
        figs <- list_mixture_plots(model2, bafdat)
        mixture_layout(figs, augmented=FALSE)    
    }

    gmodel2 <- genotype_model(model2, snp_se)
    expect_identical(modelName(model2), modelName(gmodel))
    gmodel2.obs <- gmodel2[!isSimulated(gmodel2)]
    gmodel.obs <- gmodel[!isSimulated(gmodel)]
    z_gmodel1 <- map_z(gmodel.obs)
    z_gmodel2 <- map_z(gmodel2.obs)        
    expect_true(mean(z_gmodel1 == z_gmodel2) > 0.995)

    ##additional starts helped
    ##mp <- McmcParams(burnin=200, iter=1000, nStarts=4)
    model3 <- cnv_models(mb.subsamp,
                         rowRanges(cnp_se)["CNP_029"],
                         snp_se,
                         mp)
    expect_identical(modelName(model3), modelName(model2))
    z_gmodel3 <- map_z(model3[ !isSimulated(model3) ])
    expect_true(mean(z_gmodel3 == z_gmodel2) > 0.995)
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


.test_that("cnp240", {
    ##
    ## This is a 37kb region with  78 SNPs
    ##
    ## It is likely that the CNV region is not well defined and that many of the SNPs fall outside the deletion boundaries
    ## 
    library(panc.data)
    data(cnp_se, package="panc.data")
    data(snp_se, package="panc.data")
    snps <- subsetByOverlaps(snp_se, cnp_se["CNP_240", ])
    full.data <- median_summary(snps,
                                assay_index=2,
                                provisional_batch=cnp_se$Sample.Plate,
                                THR=-1)
    batched.data <- kolmogorov_batches(full.data, 1e-6)
    set.seed(1234)
    downsampled.data <- down_sample2(batched.data, min_size=250)
    mb <- MultiBatch("MB3", data=downsampled.data)
    expect_true(validObject(mb))
    mp <- McmcParams(burnin=100, iter=1000, thin=1)
    devtools::load_all()
    fit <- homdeldup_model(mb, mp)
    ggMixture(fit)
    ggMixture(fit[ !isSimulated(fit) ], bins=200 )
    bmodel <- fit[!isSimulated(fit)]
    expect_true(all(id(bmodel) %in% colnames(snps)))
    gmodel <- genotype_model(bmodel, snps)
    tmp <- upsample2(gmodel, full.data)
    freq <- as.integer(table(tmp$copynumber))
    gap::hwe(freq, data.type="count")
    if(FALSE){
        model2 <- select_highconfidence_samples2(gmodel, snps)
        bafdat <- join_baf_oned(gmodel, snps)
        figs <- list_mixture_plots(model2, bafdat)
        mixture_layout(figs, augmented=FALSE)
    }
})
