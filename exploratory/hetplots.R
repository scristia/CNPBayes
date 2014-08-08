library(snow)
library(doSNOW)
library(foreach)
library(oligoClasses)
library(CNPBayes)
library(aricExperimentData)
library(aricUricAcid)
library(GenomicRanges)
expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"
NW <- 15

percent.het <- function(b) {
    is.het <- ifelse(b <= 0.6 & b >= 0.4, 1, 0)
    phet <- mean(is.het, na.rm=TRUE)
}

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]

vi.gr <- readRDS(file="VI-granges.rds")
names(vi.gr) <- paste0("CNP", seq_along(vi.gr))
if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(chromosome(feature.gr) %in% autosomes())
feature.gr <- feature.gr[auto.index]

snp.ind <- feature.gr$is.snp == TRUE
wc.files <- list.files(expdd, full.names=TRUE, pattern="wc_")
names(wc.files) <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
wc.files <- wc.files[as.character(pd.df$cel_id)] #8411

hits <- findOverlaps(vi.gr, feature.gr)
hitlist <- split(subjectHits(hits), queryHits(hits))
index <- subjectHits(hits)

## count number of SNPs
snpcnt <- sapply(1:length(hitlist), function(x) sum(snp.ind[hitlist[[x]]]))

batchjob <- rep(seq_len(NW), length.out=length(wc.files))
fileList <- split(wc.files, batchjob)
fileList <- fileList[1:NW]

het <- foreach(k = seq_along(fileList), .combine = rbind) %dopar% {
    files <- fileList[[k]]
    MAT <- matrix(NA, length(files), length(vi.gr))
    rownames(MAT) <- names(files)
    for(m in seq_along(files)) {
        dat <- readRDS(wc.files[m])
        b <- dat[, 3]
        b <- b/1000
        b[!snp.ind] <- NA
        blist <- sapply(1:length(hitlist), function(x) b[hitlist[[x]]])
        #blist <- split(b, queryHits(hits))

        MAT[m, ] <- sapply(blist, percent.het)
    }
    MAT
}
colnames(het) <- names(vi.gr)
## create matrix: sample by loci
saveRDS(het, file="baf-mat.rds")
q('no')
