library(CNPBayes)
library(GenomicRanges)
library(oligoClasses)
library(aricExperimentData)
library(aricUricAcid)

expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]

## EA samples
pd.ea <- pd.df[EA(pd.df), ]
pd.ea <- pd.ea[keep(pd.ea), ]

## AA samples
pd.aa <- pd.df[AA(pd.df), ]
pd.aa <- pd.aa[keep(pd.aa), ]
rm(pd.df)

### EA with quality controls
files <- list.files(expdd, pattern="vi_")
names(files) <- gsub(".rds", "", gsub("vi_", "", files))
files <- files[rownames(pd.ea)]

vi.grl <- readOrWriteFile(file.path(uadir, "CNV_grl_EA_VI.rds"),
                          collectGRanges,
                          datadir=expdd,
                          files=files,
                          max.cnv=100,
                          force=FALSE) ##8955
#penn.grl <- readOrWriteFile("/home/student/scristia/Software/CNPBayes/penncnv_granges.rds",
#                          collectGRanges,
#                          datadir=expdd,
#                          files=files,
#                          max.cnv=100,
#                          force=FALSE) ##8955

lrr.stats <- readRDS(file.path(expdd, "lrr_stats_penn.rds"))
lrr.stats <- lrr.stats[names(vi.grl), ]
drop <- which(elementLengths(vi.grl) > 100 | lrr.stats[, "lag10ACF"] > 0.03 | lrr.stats[, "mad"] > 0.32)
vi.grl <- vi.grl[-drop]
pd.ea <- pd.ea[names(vi.grl), ]
rm(lrr.stats, drop)

ea.names <- rownames(pd.ea)
lrr.expdata <- readRDS("~/Software/summarized_experiment_average_lrr_hg18.rds")
#lrr.expdata <- readRDS("~/Software/CNPBayes/se_median_lrr_penncnv.rds")
n <- length(rowData(lrr.expdata))
ea.inds <- match(ea.names, colnames(lrr.expdata))

avgRs <- sapply(1:n, function(x) assays(lrr.expdata[x])$avglrr[ea.inds]/100)
rownames(avgRs) <- colnames(lrr.expdata)[ea.inds]
saveRDS(avgRs, "~/Software/avgRs_wc_ea-vi.rds")
#saveRDS(avgRs, "~/Software/avgRs_wc_ea-penn.rds")

aa.names <- rownames(pd.aa)
aa.inds <- match(aa.names, colnames(lrr.expdata))
avgRs.aa <- sapply(1:n, function(x) assays(lrr.expdata[x])$avglrr[aa.inds]/100)
rownames(avgRs.aa) <- colnames(lrr.expdata)[aa.inds]
saveRDS(avgRs.aa, "~/Software/avgRs_wc_aa_NOQC-vi.rds")
#saveRDS(avgRs.aa, "~/Software/avgRs_wc_aa_NOQC-penn.rds")




########### mccarroll stuff
mcrl <- readRDS("se_mccarroll.rds")
n <- length(rowData(mcrl))
ea.inds <- match(ea.names, colnames(mcrl))
rs <- sapply(1:n, function(x) assays(mcrl[x])$avglrr[ea.inds]/100)
rownames(rs) <- colnames(mcrl)[ea.inds]

pdf('~/tmp/mccarrollcnp_plots%03d.pdf', onefile=FALSE)
for(i in 1:ncol(rs)) {
    main <- paste0("chr1: ", start(mcrl@rowData[i]), "-", end(mcrl@rowData[i]))
    tmp <- hist(rs[,i], breaks=200, col="lightgray", border="lightgray", freq=FALSE, main=main)
    text(tmp$breaks[1], max(tmp$density), lab=rowData(mcrl)$CN[i])
}
dev.off()
