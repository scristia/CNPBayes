set.seed(9209)
library(SummarizedExperiment)
library(panc.data)
data(cnp_se, package="panc.data")
data(snp_se, package="panc.data")
colnames(cnp_se) <- seq_along(colnames(cnp_se))
colData(cnp_se) <- DataFrame(plate=colData(cnp_se)[, "Sample.Plate"])
cnp_se <- cnp_se[c(1, 14, 22, 29, 147, 244), ]
rowRanges(cnp_se) <- granges(rowRanges(cnp_se))
medians <- assays(cnp_se)[["MEDIAN"]]
se <- SummarizedExperiment(assays=SimpleList(MEDIAN=medians),
                           rowRanges=rowRanges(cnp_se),
                           colData=colData(cnp_se))
colnames(se) <- seq_len(ncol(se))
cnp_se <- se

snp_se <- subsetByOverlaps(snp_se, cnp_se)
saveRDS(snp_se, file="../inst/extdata/snp_se.rds")
saveRDS(cnp_se, file="../inst/extdata/cnp_se.rds")
