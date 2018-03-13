grl <- readRDS(system.file("extdata", "grl_deletions.rds", package="CNPBayes"))
cnv.region <- consensusCNP(grl, max.width=5e6)
saveRDS(cnv.region, file="../inst/extdata/cnv_region.rds")
