context("Identify batches")

test_that("summarize_region", {
    library(SummarizedExperiment)
    set.seed(2463)
    ##
    ## summarize_region code chunk
    ##
    path <- system.file("extdata", package="CNPBayes")
    cnp_se <- readRDS(file.path(path, "cnp_se.rds"))
    plates <- colData(cnp_se)$Sample.Plate
    ##
    ## To keep the unit test short, focus on first 20 plates
    ##
    plates1_20 <- plates[ plates %in% unique(plates)[1:20] ]
    keep <- colData(cnp_se)$Sample.Plate %in% plates1_20
    dat <- summarize_region(cnp_se[1, keep],
                            provisional_batch=plates1_20,
                            THR=-1,
                            min_size=300)
    expect_identical(length(unique(dat$batch)), 5L)
    expect_identical(nrow(dat), 1222L)
    if(FALSE){
        ##
        ## Save downsampled data from entire dataset
        ##
        set.seed(14)
        dat <- summarize_region(cnp_se[1, ],
                                provisional_batch=cnp_se$Sample.Plate,
                                THR=-1,
                                S=1500,
                                min_size=200,
                                KS_cutoff=1e-5)
        group_by(dat, batch) %>%
            tally()

        ## If S is the size of the dataset, their is no down-sampling
        set.seed(14)        
        dat2 <- summarize_region(cnp_se[1, ],
                                 provisional_batch=cnp_se$Sample.Plate,
                                 THR=-1,
                                 S=6038,
                                 min_size=200,
                                 KS_cutoff=1e-5)
        tal <- group_by(dat2, batch) %>%
            tally()
        expect_identical(sum(tal$n), ncol(cnp_se))
        saveRDS(dat, file=file.path("..",
                                    "..",
                                    "inst",
                                    "extdata",
                                    "CNP_001",
                                    "batched_data.rds"))
        library(ggplot2)
        ggplot(dat, aes(oned)) +
            geom_histogram(bins=200) +
            facet_wrap(~batch, ncol=1)
    }
})
