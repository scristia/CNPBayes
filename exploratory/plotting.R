source('segplots.R')
library(CNPBayes)
library(oligoClasses)
library(GenomicRanges)
library(aricExperimentData)
library(aricUricAcid)

expdd <- "/local_data/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"

reduced.gr <- readRDS('~/Software/wcreduced.gr')
avgRs <- readRDS('~/Software/avgRs_wc_ea.rds')
vi.grl <- readRDS(file.path(uadir, "CNV_grl_EA_VI.rds"))
gr <- unlist(vi.grl)

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(chromosome(feature.gr) %in% autosomes())
feature.gr <- feature.gr[auto.index]

olaps <- findOverlaps(reduced.gr, gr)
indices <- split(subjectHits(olaps), queryHits(olaps))

## find marker positions
ix <- findOverlaps(feature.gr, reduced.gr)
mx <- split(queryHits(ix), subjectHits(ix))

## LRR list
#lrr.stats <- readRDS(file.path(uadir, "lrr_stats_penn.rds"))
#lrr.list <- lapply(1:length(mx), function(x) feature.gr[mx[[x]]])

## mixture model params
delta=0.1
S=1000
burnin=200

## init matrix for diffs/variances
matx <- matrix(NA, nrow=length(avgRs[1,]), ncol=2)
pdf('~/tmp/wcplots%03d.pdf', width=8, height=8, onefile=FALSE)
for(i in seq_along(avgRs[1,])) {
    markers <- feature.gr[mx[[i]]]
    par(mar = c(0, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
    layout(matrix(c(1,3,2,4), 2, 2, byrow=TRUE))
    segplots(indices=indices[[i]], red.range=reduced.gr[i], subject.cnvs=gr,
             flag=100000L, olap=TRUE, markers=markers)
    posts <- getPost(r=avgRs[,i], kmin=1, kmax=5, delta=delta, S=S,
                     burnin=burnin, plot=TRUE)

    hist(avgRs[,i], breaks=200, col="lightgray", border="lightgray",
         freq=FALSE, main="")
    K <- which.min(sapply(1:5, function(x) posts[[x]]$bic))
    PI <- colMeans(posts[[K]]$P)
    MU <- colMeans(posts[[K]]$means)
    SIGMA <- colMeans(1/posts[[K]]$precs)
    y <- seq(min(avgRs[,i]), max(avgRs[,i]), len=1000)
    for(k in 1:K) lines(y, PI[k]*dnorm(y, MU[k], sqrt(SIGMA[k])),
                        col='gray40', lwd=2)
    ## name plot
    region <- paste0("Region ",i)
    chr <- chromosome(reduced.gr[i])
    mtext(text=paste(chr, ": ", region), side=3, outer=TRUE)
    if(K > 1) {
        matx[i,1] <- diff(MU)[1]
        matx[i,2] <- SIGMA[1]
    }
}
dev.off()
saveRDS(matx, "diffvar_matrix.RDS")
q('no')
