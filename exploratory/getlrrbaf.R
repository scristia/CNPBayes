library(snow)
library(doSNOW)
library(foreach)
library(GenomicRanges)
library(oligoClasses)
library(aricExperimentData)
library(aricUricAcid)

NW <- 20
setParallelization(WORKERS=NW)

expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
#pd.df <- pd.df[EA(pd.df), ]
pd.df <- pd.df[keep(pd.df), ]

wc.files <- list.files(expdd, full.names=TRUE, pattern="wc_")
names(wc.files) <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
wc.files <- wc.files[as.character(pd.df$cel_id)] #8411

#wc.files <- sample(wc.files, 20) ## for testing

batchjob <- rep(seq_len(NW), length.out=length(wc.files))
fileList <- split(wc.files, batchjob)

### which region to consider
chr5.gr <- GRanges("chr5", IRanges(8755522, 8800142))

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
chr5.index <- which(chromosome(feature.gr) == "chr5")
feature.gr <- feature.gr[chr5.index]


olaps <- findOverlaps(chr5.gr, feature.gr)
subj.index <- subjectHits(olaps)
query.index <- queryHits(olaps)

fileList <- fileList[1:NW]
for(i in seq_len(NW)) fileList[[i]] <- fileList[[i]]
lrrbaf <- foreach(k = seq_along(fileList), .combine="cbind") %dopar% {
    files <- fileList[[k]]
    A <- matrix(NA, ncol=length(files))
    colnames(A) <- gsub(".rds", "", gsub("wc_", "", basename(files)))
    for(m in seq_along(files)){
        dat <- readRDS(files[m])
        rlist <- split(dat[subj.index, "robust_lrr"], query.index)
#        baflist <- split(dat[subj.index, "BAF"], query.index)
        A[,m] <- rlist[[1]]
    }
    A
}



## check low level data
xfile <- "/local/recover/r00/aric/aricExperimentData/data/wc_243080.rds"
xdat <- readRDS(xfile)
marker.index <- subjectHits(findOverlaps(chr5.gr, feature.gr, maxgap=500e3))
#atv.lrr <- xdat[marker.index, "atv_lrr"]/100
robust.lrr <- xdat[marker.index, "robust_lrr"]/100
baf <- xdat[marker.index, "BAF"]/1000
baf[baf==2] <- 0

pos <- start(feature.gr[marker.index])

xlim <- range(pos)
par(mfrow=c(2,1), las=1, mar=c(0.1,4,0.1,0.1), oma=c(4,0.1,0.1,0.1))
#plot(pos, atv.lrr, pch=20, xaxt="n", col="gray50")
#abline(v=c(start(chr5.gr), end(chr5.gr)), lty=2)
plot(pos, robust.lrr, pch=20, xaxt="n", col="gray50")
abline(v=c(start(chr5.gr), end(chr5.gr)), lty=2)
plot(pos, baf, pch=20, xaxt="n", col="gray50")
abline(v=c(start(chr5.gr), end(chr5.gr)), lty=2)
at <- pretty(xlim, n=8)
axis(1, at=at, labels=at/1e6, cex.axis=0.7)

