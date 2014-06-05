## misc testing stuff
#avgRs <- readRDS('~/Software/avgRs_wc_ea.rds')
avgRs <- readRDS('../CNPBayes/data/avgRs_wc_ea.rds')
rquant <- quantile(r, c(0.01, 0.99))
#rr <- r[r > rquant[1] & r < rquant[2]]
findpeaks <- function(r, delta, main="") {
    rquant <- quantile(r, c(0.01, 0.99))
    rr <- r[r > rquant[1] & r < rquant[2]]
    rd <- density(rr, adjust=1)
    ## find inflection points (sign changes from + to -)
    itx <- which((diff(sign(diff(rd$y))) == -2)) + 1
    ## 50th quantile of y is very ad hoc but seems to work well
#    ity <- which(rd$y[itx] > quantile(rd$y, 0.5))
#    ity <- which(rd$y[itx] > 0.1)
    ity <- which(rd$y[itx] > 0.01*max(rd$y))
    itx <- itx[ity]

    ## short algorithm for combining nearby peaks
    modes <- rd$x[itx]
    yval <- rd$y[itx]
    v <- as.integer(diff(modes) < delta)
    v <- append(v, 0L)
    peaks <- rep(NA, sum(v==0))
    vrev <- rev(v)
    modesrev <- rev(modes)
    yvalr <- rev(yval)
    zind <- which(vrev == 0)
    for(i in seq_along(zind)) {
#        if(i < length(peaks)) peaks[i] <- mean(modesrev[zind[i]:(zind[i+1]-1L)])
#        else peaks[i] <- mean(modesrev[zind[i]:length(modesrev)])
        if(i < length(peaks)) {
            ymax <- which.max(yvalr[zind[i]:(zind[i+1]-1L)])
            peaks[i] <- modesrev[zind[i]:(zind[i+1]-1L)][ymax]
        }
        else {
            ymax <- which.max(yvalr[zind[i]:length(yvalr)])
            peaks[i] <- modesrev[zind[i]:length(modesrev)][ymax]
        }
    }
    peaks <- rev(peaks)

    kmpeaks <- kmeans(rr, length(peaks), nstart=10)$centers
    hist(r, breaks=200, col=rgb(211,211,211,alpha=150,max=255), border=rgb(211,211,211, max=255), freq=FALSE, main=main)
    lines(density(r), lwd='2', col=rgb(0, 153, 204, max=255))
    abline(v=peaks, col='gray40', lwd='2')
    abline(v=kmpeaks, col='skyblue3', lwd='2', lty=2)
    abline(v=modes, lty=2)
#    abline(v=modes[-peaks], col='gray40', lwd='2')
#   abline(h=quantile(rd$y,.5), col="orange")
}

pdf('~/scratch/plots/findpeaks%03d.pdf', onefile=FALSE)
par(mfrow=c(2,2))
for(j in seq_along(avgRs[155,])) {
    findpeaks(r = avgRs[,j], delta = 0.15, main=paste("Region", j))
}
dev.off()


### test
r <- avgRs[,155]
S = 1000
delta = 0.15
burnin=0
posts <- getPost(r1, kmin=1, kmax=5, delta=delta, S=S, burnin=burnin, plot=TRUE)
