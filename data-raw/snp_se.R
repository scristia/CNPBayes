library(panc.data)
library(CNPBayes)
##
## summarize_region code chunk
##
data(cnp_se, package="panc.data")
data(snp_se, package="panc.data")
##g <- GRanges("chr1", IRanges(1627805, 1673809),
##             seqinfo=Seqinfo("chr1", 249250621,
##                             genome="hg19"))
##cnp_se <- cnp_se[1, ]
cnp_ids <- paste0("CNP_", c("001", "014", "022", "029", "147", "244"))
cnp_se <- cnp_se[cnp_ids, ]
id <- colnames(cnp_se)
snp_se <- snp_se[overlapsAny(snp_se, cnp_se), ]
##coldat <- DataFrame(Sample.Plate=colData(cnp_se)[, "Sample.Plate"])
##colData(cnp_se) <- coldat
colnames(cnp_se) <- seq_along(id)
colnames(snp_se) <- colnames(cnp_se)
colData(snp_se) <- colData(cnp_se)
saveRDS(cnp_se, file="../inst/extdata/cnp_se.rds")
saveRDS(snp_se, file="../inst/extdata/snp_se.rds")
