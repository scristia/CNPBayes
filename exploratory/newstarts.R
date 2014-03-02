## misc testing stuff
avgRs <- readRDS('~/Software/avgRs_wc_ea.rds')
rquant <- quantile(r, c(0.01, 0.99))
#rr <- r[r > rquant[1] & r < rquant[2]]
findpeaks <- function(r, delta, main) {
    rd <- density(r, adjust=1)
    ## find inflection points (sign changes from + to -)
    itx <- which((diff(sign(diff(rd$y))) == -2)) + 1
    ## 50th quantile of y is very ad hoc but seems to work well
    ity <- which(rd$y[itx] > quantile(rd$y, 0.5))
#    ity <- which(rd$y[itx] > 0.1)
    #ity <- which(rd$y[itx] > 0.05*max(rd$y))
    itx <- itx[ity]

    ## short algorithm for combining nearby peaks (make use of RLE)
    if(any(diff(rd$x[itx]) < delta)) {
        v <-diff(rd$x[itx]) < delta
        v.rle <- rle(v)
        peaks <- rep(NA, length(rd$x[itx]))
        inc <- 1
        for(i in seq_along(v.rle$lengths)) {
            if( v.rle$values[i] == 1) {
                i1 <- sum(v.rle$lengths[1:i-1]) + 1
                i2 <- sum(v.rle$lengths[1:i]) + 1
                m <- rd$x[itx][i1:i2][which.max(rd$y[itx][i1:i2])]
#                m <- (rd$x[itx][i1] + rd$x[itx][i2])/2
                peaks[inc] <- m
                inc <- inc + 1
            }
            ## find isolated peaks
            else if(v.rle$values[i] == 0) {
                ## don't save if a gap between unique modes
                i3 <- sum(v.rle$lengths[1:i]) - v.rle$lengths[i] + 1
                for(j in 1:v.rle$lengths[i]) {
                    peaks[inc] <- rd$x[itx][i3 + j]
                    inc <- inc + 1
                }
#                if(i < length(v.rle$lengths) &
#                   diff(c(v.rle$lengths[i], v.rle$lengths[i+1])) < 0.1) {
#                    i3 <- sum(v.rle$lengths[1:i]) - v.rle$lengths[i] + 1
#                    for(j in 1:v.rle$lengths[i]) {
#                        peaks[inc] <- rd$x[itx][i3 + j - 1]
#                        inc <- inc + 1
#                    }
#                }
                if(i == 1) {
                    peaks[inc] <- rd$x[itx][1]
                    inc <- inc + 1
                }
            }
        }
        peaks <- peaks[!is.na(peaks)]
    } else peaks <- rd$x[itx]

    hist(r, breaks=200, col=rgb(211,211,211,alpha=150,max=255), border=rgb(211,211,211, max=255), freq=FALSE, main=main)
    lines(density(r), lwd='2', col=rgb(0, 153, 204, max=255))
    abline(v=peaks, col='gray40', lwd='2')
    abline(v=rd$x[itx], col='gray40', lwd='2')
#    abline(h=quantile(rd$y,.5), col="orange")
}

pdf('~/tmp/findpeaks%03d.pdf', onefile=FALSE)
par(mfrow=c(2,2))
for(j in seq_along(avgRs[155,])) {
    findpeaks(r = avgRs[,155], delta = 0.1, main=paste("Region", j))
}
dev.off()


### test
r <- avgRs[,155]
S = 1000
delta = 0.15
burnin=0
posts <- getPost(r1, kmin=1, kmax=5, delta=delta, S=S, burnin=burnin, plot=TRUE)
