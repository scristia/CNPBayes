library(CNPBayes)
avgRs <- readRDS("../../avgRs_wc_ea-vi.rds")
avgRs.aa <- readRDS("../../avgRs_wc_aa_NOQC-vi.rds")

## region 122 - check
x <- avgRs[,347]
x.aa <- avgRs.aa[,113]
#x.aa <- avgRs.aa[,186]
xx <- c(x, x.aa)

S=1500
burnin=500
delta=1e-5
#quant <- quantile(x, c(0.001, 0.999))
pdf("/Users/scrist/Documents/Presentations/cnpbayes/reg132post.pdf")
### Region 347 - bad variance
set.seed(4321)
post <- getPost(r=x, nu0=1, kappa0=1, kmin=1, kmax=4, delta=delta, S=S, tau20=0.1,
                 burnin=burnin, plot=TRUE, crit="bic", full=FALSE)
dev.off()

source('segplots.R')


assignments <- matrix(NA, nrow=length(xx), ncol=5, dimnames=list(names(xx), seq(1:5)))

avgRs.all <- rbind(avgRs, avgRs.aa)
posteriors <- list()
for(i in 1:5) {
    xx <- avgRs.all[,i]
    post <- getPost(r=xx, kmin=1, kmax=5, delta=delta, S=S, tau20=0.1,
                    burnin=burnin, plot=TRUE, crit="icl", full=FALSE)
    if(post$K > 1) {
        ind <- match(rownames(post$Z), names(xx))
        assignments[ind, i] <- apply(post$Z, 1, which.max)
    }
    K <- post$k
    PI <- colMeans(post$P)
    MU <- colMeans(post$means)
    SIGMA <- colMeans(1/post$precs)
    Z <- post$Z
#    snp.ids <- grep("SNP_", markers$feature.id)
#    snp.pos <- start(markers[snp.ids])
    posteriors <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z)
#    posts.EA[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, snps=snp.pos)

}

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
## check low level diiiiasdflkasdjfiivvbvqqwer:q
gtgt:tabe:Tabe
asdfdsfdfasdfasdfgawefaweta
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

avgRs <- readRDS("./home/student/scristia/Software/avgRs_wc_ea-vi.rds")
posts <- readRDS("/home/student/scristia/posts_vi2.rds")
xx <- avgRs[,2]

plotPosts(xx, posts[[2]])
mus <- colMeans(posts[[2]]$means)
sds<- 1/sqrt(colMeans(posts[[2]]$precs))
ps <- colMeans(posts[[2]]$P)
mix <- mixture(mus, sds, pis)
hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
y <- seq(-3, 2, len=1000)
for(k in 1:length(mix)) {
    if(k < length(mix)) {
        if(length(mix[[k]]) >= 2 & length(mix[[k+1]]) > length(mix[[k]]))
            next
    }
    d <- rowSums(sapply(mix[[k]], mix.dens, y=y,
                        p=ps, mix.mean=mus, mix.sd=sds), na.rm=TRUE)
    lines(y, d, col="gray40", lwd=2)
}

xs <- seq(min(mus[3] - 3*sds[3], mus[4] - 3*sds[4]), max(mus[3] + 3*sds[3],
                                                         mus[4] + 3*sds[4]), .01)
f1 <- ps[3]*dnorm(xs, mean=mus[3], sd=sds[3])
f2 <- ps[4]*dnorm(xs, mean=mus[4], sd=sds[4])

ys <- min.f1f2(xs, MU1=mus[3], MU2=mus[4], SD1=sds[3], SD2=sds[4], PI1=ps[3], PI2=ps[4])
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])

hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
y <- seq(min(xx), max(xx), len=1000)
for(k in 1:K) lines(y, ps[k]*dnorm(y, mus[k], sds[k]), col='gray40', lwd=2)
polygon(xs, ys, col="gray40", density=50, angle=45)
