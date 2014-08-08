source('segplots.R')
library(CNPBayes)
library(oligoClasses)
library(GenomicRanges)
library(aricExperimentData)
library(aricUricAcid)

expdd <- "/local_data/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"

## penn and vi - comment one out
reduced.gr <- readRDS('~/Software/wcreduced.gr')
#reduced.gr <- readRDS('~/Software/wcreduced.gr')
#reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')
## EA
avgRs <- readRDS('~/Software/avgRs_wc_ea-vi.rds')
#avgRs <- readRDS('~/Software/avgRs_wc_ea-penn.rds')
## AA
avgRs.aa <- readRDS('~/Software/avgRs_wc_aa_NOQC-vi.rds')
#avgRs.aa <- readRDS('~/Software/avgRs_wc_aa_NOQC-penn.rds')
vi.grl <- readRDS(file.path(uadir, "CNV_grl_EA_VI.rds"))
#penn.grl <- readRDS('~/Software/CNPBayes/penncnv_granges.rds')
gr <- unlist(vi.grl)
#gr <- unlist(penn.grl)

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(chromosome(feature.gr) %in% autosomes())
feature.gr <- feature.gr[auto.index]

olaps <- findOverlaps(reduced.gr, gr)
indices <- split(subjectHits(olaps), queryHits(olaps))

## find marker positions
ix <- findOverlaps(feature.gr, reduced.gr, maxgap=1e6)
mx <- split(queryHits(ix), subjectHits(ix))

## LRR list
#lrr.stats <- readRDS(file.path(uadir, "lrr_stats_penn.rds"))
#lrr.list <- lapply(1:length(mx), function(x) feature.gr[mx[[x]]])

## mixture model params
delta=0.15
S=1500
burnin=500


avgRs.all <- rbind(avgRs, avgRs.aa)
assignments <- matrix(NA, nrow=nrow(avgRs.all), ncol=ncol(avgRs.all),
                      dimnames=list(rownames(avgRs.all), seq(1:ncol(avgRs.all))))
posteriors <- list()
length(posteriors) <- length(avgRs.all[1,])

pdf('~/tmp/vicnp-plots%03d.pdf', width=11, height=8, onefile=FALSE)
for(i in seq_along(avgRs[1,])) {
    print(paste('REGION',i))
    markers <- feature.gr[mx[[i]]]
    m2 <- markers[queryHits(findOverlaps(markers, reduced.gr[i]))]
    par(mar = c(0, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
#    layout(matrix(c(1,3,5,2,4,6), 2, 3, byrow=TRUE))
    layout(matrix(c(1,2), 1, 2, byrow=TRUE))
    segplots(indices=indices[[i]][-183], red.range=reduced.gr[i], subject.cnvs=gr,
             flag=100000L, olap=FALSE, markers=markers)

    ### EUROPEAN ANCESTRY
    xx <- avgRs.all[,i]
    set.seed(999)
    posts <- getPost(r=xx, kmin=1, kmax=5, delta=delta, S=S,
                     burnin=burnin, plot=TRUE, crit="icl")

    K <- posts$K
    PI <- colMeans(posts$P)
    MU <- colMeans(posts$means)
    SIGMA <- colMeans(1/posts$precs)
    Z <- posts$Z
    snp.ids <- grep("SNP_", m2$feature.id)
    snp.pos <- start(m2[snp.ids])
    if(K > 1) {
        ind <- match(rownames(posts$Z), names(xx))
        assignments[ind, i] <- apply(posts$Z, 1, which.max)
    }
    posteriors[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, snps=snp.pos)
#    posts.EA[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, snps=snp.pos)
#    y <- seq(min(avgRs[,i]), max(avgRs[,i]), len=1000)
#    for(k in 1:K) lines(y, PI[k]*dnorm(y, MU[k], sqrt(SIGMA[k])),
#                        col='gray40', lwd=2)


    ## AFRICAN AMERICAN ANCESTRY
#    posts <- getPost(r=avgRs.aa[,i], kmin=1, kmax=5, delta=delta, S=S,
#                     burnin=burnin, plot=TRUE, crit="icl")
#
##    hist(avgRs.aa[,i], breaks=200, col="lightgray", border="lightgray",
##         freq=FALSE, main="")
#    K <- which.min(sapply(1:5, function(x) posts[[x]]$icl))
#    PI <- colMeans(posts[[K]]$P)
#    MU <- colMeans(posts[[K]]$means)
#    SIGMA <- colMeans(1/posts[[K]]$precs)
#    Z <- posts[[K]]$Z
#    snp.ids <- grep("SNP_", markers$feature.id)
#    snp.pos <- start(markers[snp.ids])
#    posts.AA[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, snps=snp.pos)
    ## name plot
    region <- paste0("Region ",i)
    chr <- chromosome(reduced.gr[i])
    mtext(text=paste(chr, ": ", region), side=3, outer=TRUE)
}

dev.off()
saveRDS(posteriors, "posteriors_all-vi.rds")
saveRDS(assignments, "Assignment-matrix.rds")
q('no')


#set.seed(999)
#posts <- invisible(getPost(r=avgRs[,155], kmin=1, kmax=5, delta=delta, S=S, burnin=0,
#        plot=TRUE, crit="icl"))
