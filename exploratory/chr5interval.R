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
avgLrr <- foreach(k = seq_along(fileList), .combine="append") %dopar% {
    files <- fileList[[k]]
    avgRs <- rep(NA, length(files))
    names(avgRs) <- gsub(".rds", "", gsub("wc_", "", basename(files)))
    for(m in seq_along(files)){
        dat <- readRDS(files[m])
        rlist <- split(dat[subj.index, "robust_lrr"], query.index)
        avgRs[m] <- sapply(rlist, median, na.rm=TRUE)/100
    }
    avgRs
}

avgRs.ea <- avgLrr[match(EA(pd.df), names(avgLrr))]
saveRDS(avgRs.ea, file="avgrs_ch5_EA.rds")

avgRs.aa <- avgLrr[match(AA(pd.df), names(avgLrr))]
saveRDS(avgRs.aa, file="avgrs_ch5_AA.rds")

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

library(CNPBayes)
## mixture model params
delta=0.15
S=1500
burnin=500

### EUROPEAN ANCESTRY
avgRs.adj <- readRDS("~/MixtureModels/data/avgRs_adj.rds")
posts <- getPost(r=avgRs.adj[,128], kmin=1, kmax=5, delta=delta, S=S,
                 burnin=burnin, plot=TRUE, crit="icl")

#    hist(avgRs[,i], breaks=200, col="lightgray", border="lightgray",
#         freq=FALSE, main="")
K <- which.min(sapply(1:5, function(x) posts[[x]]$icl))
PI <- colMeans(posts[[K]]$P)
MU <- colMeans(posts[[K]]$means)
SIGMA <- colMeans(1/posts[[K]]$precs)
Z <- posts[[K]]$Z
assignments <- apply(Z, 1, which.max)

quant <- quantile(avgRs.adj[,128], c(0.001, 0.999))
samples <- names(avgRs.adj[,128][avgRs.adj[,128] > quant[1] & avgRs.adj[,128] < quant[2] ])
posts.EA <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, assignments=assignments,
                 samples=samples)
saveRDS(posts.EA, "chr5reg_EA.rds")

write.table(cbind(as.integer(samples), assignments), "chr5_assignments.txt", row.names=FALSE, col.names=FALSE)

q('no')
