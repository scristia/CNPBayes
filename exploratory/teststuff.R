library(CNPBayes)
avgRs <- readRDS("../avgRs_wc_ea-vi.rds")
avgRs.aa <- readRDS("../avgRs_wc_aa_NOQC-vi.rds")
## region 122 - check
x <- avgRs[,140]
x.aa <- avgRs.aa[,140]
#x.aa <- avgRs.aa[,186]
xx <- c(x, x.aa)

S=1500
burnin=500
delta=0.1
#quant <- quantile(x, c(0.001, 0.999))
post <- getPost(r=xx, kmin=4, kmax=5, delta=delta, S=S, tau20=0.1,
                 burnin=burnin, plot=TRUE, crit="icl", full=FALSE)


### AFRICAN AMERICANS
set.seed(219)
posts.aa <- getPost(r=x.aa, kmin=4, kmax=5, delta=delta, S=S, tau20=0.1,
                 burnin=burnin, plot=TRUE, crit="icl", full=FALSE)


#### segfault at region 140 in african americans

#### check region 186
post186 <- readRDS("posteriors_EA.RDS")
post186 <- post186[[186]]

post186.aa <- readRDS("posteriors_AA.RDS")
post186.aa <- post186.aa[[186]]

########################################################################
library(CNPBayes)
library(oligoClasses)
library(GenomicRanges)
library(aricExperimentData)
library(aricUricAcid)

expdd <- "/local_data/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"


#reduced.gr <- readRDS('~/Software/wcreduced.gr')
reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')

#vi.grl <- readRDS(file.path(uadir, "CNV_grl_EA_VI.rds"))
penn.grl <- readRDS('~/Software/CNPBayes/penncnv_granges.rds')
#gr <- unlist(vi.grl)
gr <- unlist(penn.grl)



library(GenomicRanges)
## check low level data
loci <- which(start(gr) == 105820728)
grx <- gr[loci]
rs <- avgRs[,paste(loci)]
delta=0.15
S=1000
burnin=200
tau20 <- 0.1
posts <- getPost(r=rs, kmin=2, kmax=3, delta=delta, S=S,
                 burnin=burnin, plot=TRUE, crit="icl", main="")
quant <- quantile(rs, c(0.001, 0.999))
rs=rs[rs > quant[1] & rs < quant[2] ]
vvv <- posts[[2]]$Z
rownames(vvv) <- names(rs)
hdels <- names(vvv[,1][vvv[,1] > 900])
ssmp <- sample(hdels, 10)

for(j in 1:10) {
    name <- paste0("~/tmp/sample",ssmp[j],".pdf")
    pdf(name, onefile=FALSE)

    xfile <- paste0("/local/recover/r00/aric/aricExperimentData/data/wc_",ssmp[j],".rds")
    xdat <- readRDS(xfile)
    marker.index <- subjectHits(findOverlaps(grx, feature.gr, maxgap=50e3))
    robust.lrr <- xdat[marker.index, "robust_lrr"]/100
    baf <- xdat[marker.index, "BAF"]/1000
    baf[baf==2] <- 0

    pos <- start(feature.gr[marker.index])

    xlim <- range(pos)
    par(mfrow=c(2,1), las=1, mar=c(0.1,4,0.1,0.1), oma=c(4,0.1,0.1,0.1))
    plot(pos, robust.lrr, pch=20, xaxt="n", col="gray50")
    abline(v=c(start(grx), end(grx)), lty=2)
    plot(pos, baf, pch=20, xaxt="n", col="gray50")
    abline(v=c(start(grx), end(grx)), lty=2)
    at <- pretty(xlim, n=8)
    axis(1, at=at, labels=at/1e6, cex.axis=0.7)
    dev.off()
}


x <- readRDS("/local/recover/r00/aric/aricExperimentData/data/234380_A.rds")
