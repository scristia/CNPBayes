library("CNPBayes")
## VI
posts <- readRDS("posts-v2.rds")
avgRs <- readRDS("../../avgRs_wc_ea-vi.rds")
## McCarroll
#posts <- readRDS("posts-mc.rds")
#avgRs <- readRDS("mc_avgRs_EA.rds")

dens.comp <- function(y, comp, p, mix.mean, mix.sd) {
    p[comp]*dnorm(y, mean=mix.mean[comp], sd=mix.sd[comp])
}
## hack of code found on stackexchange, which I can not find again
min.f1f2 <- function(x, MU1, MU2, SD1, SD2, PI1, PI2) {
    f1 <- rowSums(sapply(1:length(PI1), dens.comp, y=x,
                          p=PI1, mix.mean=MU1, mix.sd=SD1), na.rm=TRUE)
    f2 <- PI2*dnorm(x, mean=MU2, sd=SD2)
    pmin(f1, f2)
}


mergecomponent <- function(x, MU, SD, PI) {
    d <- diff(MU)/sqrt(sum(SD^2))
    area <- integrate(min.f1f2, -Inf, Inf, MU1=MU[1], MU2=MU[2], SD1=SD[1], SD2=SD[2], PI1=PI[1], PI2=PI[2])$value
    v <- max(area/PI)
    return(c(d, v))
}

mixture <- function(MU, SD, PI) {
    mixtures <- vector("list", K)
    i1 <- 1
    flag <- TRUE
    for(i in 1:(K-1)) {
        area <- integrate(min.f1f2, -Inf, Inf, MU1=MU[i1:i], MU2=MU[(i+1)],
                          SD1=SD[i1:i], SD2=SD[(i+1)], PI1=PI[i1:i], PI2=PI[(i+1)],
                          subdivisions=500L)$value
        v <- max(area/c(sum(PI[i1:i]), PI[(i+1)]))

        ### if area is greater than 0.5, treat as mixture
        if(v >= 0.5) {
            mixtures[[i]] <- i1:(i+1)
            if(length(i1:(i+1)) >= 3) mixtures[[i-1]] <- NULL
            flag <- TRUE
        }

        else if(v < 0.5) {
            if(flag == FALSE & i < K-1) {
                mixtures[[i]] <- i
                   }
            if(i < K-1 & flag == TRUE) {
                if(i == 1) mixtures[i] <- i
                mixtures[[i+1]] <-  i+1
                flag = FALSE
            }
            i1 <- i+1
            if(i == K-1) {
                if(flag == FALSE | i == 1) mixtures[[i]] <-  i
                mixtures[[i+1]] <-  i1
            }
        }
    }
    mixtures <- mixtures[!sapply(mixtures, is.null)]
    return(mixtures)
}

#m <- 114
#pdf('~/tmp/merge-plots%03d.pdf', width=12, height=8, onefile=FALSE)
pdf('~/tmp/dpp-mc%03d.pdf', width=12, height=8, onefile=FALSE)
mixs <- vector("list", length=length(posts))
for(m in 1:ncol(avgRs)) {
    xx <- avgRs[,m]
    MU <- colMeans(posts[[m]]$means[500:1500,])
    ind <- order(MU)
    SD <- 1/sqrt(colMeans(posts[[m]]$precs[500:1500,]))
    PI <- colMeans(posts[[m]]$P[500:1500,])
    MU <- MU[ind]
    SD <- SD[ind]
    PI <- PI[ind]
    K <- 5
#    K <- posts[[m]]$K
#    if(length(MU) == 1) {
#        mixs[[m]] <- 1
#        y <- seq(-3, 2, len=1000)
#        par(mar = c(1, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
#        layout(matrix(c(1,2), 1, 2, byrow=TRUE))
#        hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
#        for(k in 1:K) lines(y, PI[k]*dnorm(y, MU[k], SD[k]),
#                            col='gray40', lwd=2)
#        hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
#        for(k in 1:K) lines(y, PI[k]*dnorm(y, MU[k], SD[k]),
#                            col='gray40', lwd=2)
#        next
#    }

    mixtures <- mixture(MU, SD, PI)

    mixs[[m]] <- mixtures

    diffs <- c()
    areas <- c()
    for(i in 1:(K-1)) {
        b <- c(mergecomponent(xx, MU[i:(i+1)], SD[i:(i+1)], PI[i:(i+1)]))
        diffs <- c(diffs, b[1])
        areas <- c(areas, b[2])
    }
    y <- seq(-3, 2, len=1000)
    par(mar = c(1, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
    layout(matrix(c(1,2), 1, 2, byrow=TRUE))
    hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
    lines(y, rowSums(sapply(1:5, function(x) PI[x]*dnorm(y, MU[x], SD[x]))))
#    plotPosts(xx, posts[[m]], main=paste("Region", m))
    legend('topleft', legend =paste(seq(1,K-1),"=", round(diffs,3), ", ", round(areas, 3), collapse="\n"), bty="n")

    hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
    for(k in 1:length(mixtures)) {
        if(k < length(mixtures)) {
            if(length(mixtures[[k]]) >= 2 & length(mixtures[[k+1]]) > length(mixtures[[k]]))
                next
        }
        d <- rowSums(sapply(mixtures[[k]], dens.comp, y=y,
                            p=PI, mix.mean=MU, mix.sd=SD), na.rm=TRUE)
        lines(y, d, col="gray40", lwd=2)
    }
}
dev.off()

### to merge the components in the posterior output
for(j in 1:length(posts)) {
    k <- length(mixs[[j]])
    posts[[j]]$Z2 <- sapply(1:k, function(x) rowSums(posts[[j]]$Z[,mixs[[j]][[x]], drop=FALSE]))
}

saveRDS(posts, "posts_vi2.rds")
#doesn't work how I want
#posts <- mapply(list, posts, "mix"=mixs, SIMPLIFY=FALSE)

########
mergecomponent(xx, MU[i:(i+1)], SD[i:(i+1)], PI[i:(i+1)])

xs <- seq(min(MU[1] - 3*SD[1], MU[2] - 3*SD[2]), max(MU[1] + 3*SD[1], MU[2] + 3*SD[2]), .01)
f1 <- PI[1]*dnorm(xs, mean=MU[1], sd=SD[1])
f2 <- PI[2]*dnorm(xs, mean=MU[2], sd=SD[2])

ys <- min.f1f2(xs, mu1=MU[1], mu2=MU[2], sd1=SD[1], sd2=SD[2], pi1=PI[1], pi2=PI[2])
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])

y <- seq(min(xx), max(xx), len=1000)
for(k in 1:K) lines(y, PI[k]*dnorm(y, MU[k], SD[k]),
                    col='gray40', lwd=2)
polygon(xs, ys, col="gray40", density=50, angle=45)

m### this works in general
integrate(min.f1f2, -Inf, Inf, mu1=MU[3], mu2=MU[4], sd1=SD[3], sd2=SD[4], pi1=PI[3], pi2=PI[4])

##############
m <- 114
xx <- avgRs[,m]
MU <- colMeans(posts[[m]]$means)
SD <- 1/sqrt(colMeans(posts[[m]]$precs))
PI <- colMeans(posts[[m]]$P)
K <- posts[[m]]$K
hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE)
lines(y, PI[1]*dnorm(y, MU[1], SD[1]) + PI[2]*dnorm(y, MU[2], SD[2]))
lines(y, PI[3]*dnorm(y, MU[3], SD[3]) + PI[4]*dnorm(y, MU[4], SD[4]))
lines(y, PI[5]*dnorm(y, MU[5], SD[5]))



#####
min.f1f2 <- function(x, mu1, mu2, sd1, sd2, pi1, pi2) {
    f1 <- pi1*dnorm(x, mean=mu1, sd=sd1)
    f2 <- pi2*dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
}

##### filter VI regions out of MC regions
library(GenomicRanges)
mc.gr <- readRDS("mc-granges.rds")
reduced.gr <- readRDS("/home/student/scristia/Software/wcreduced.gr")
mc.gr2 <- mc.gr[match(colnames(avgRs), names(mc.gr))]

hits <- subjectHits(findOverlaps(reduced.gr, mc.gr2))
mc.gr3 <- mc.gr2[-hits]
posts2 <- posts[-hits]
saveRDS(posts2, "posts_mc-filtered.rds")
saveRDS(mc.gr3, "granges_mc-filtered.rds")

### test DP stuff
xx <- avgRs[,91]
posts <- readRDS("postdpp.rds")
dp <- dppgibbs(xx, H=5, alpha=1, S=1500)
MU <- colMeans(dp$means[500:1500,])
SD <- 1/sqrt(colMeans(dp$precs[500:1500,]))
PI <- colMeans(dp$P[500:1500,])
K <- 5
hist(xx, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main="",  xlim=c(-3, 2))
lines(y, rowSums(sapply(1:5, function(x) PI[x]*dnorm(y, MU[x], SD[x]))))
