## Make function that creates a plot for a single region
## What is passed into function:
## -reduced, indices, avgRs (r), plates, avgRs, avgRs
plotRegion <- function(r, r.adj, indices, red.range, adj.range, subject.cnvs,
                       df, region, markers, markers.adj) {
    j <- indices
    x <- start(subject.cnvs)[j]
    y <- end(subject.cnvs)[j]
    ranges <- cbind(x, y)

    w <- cbind(red.range[1] - ranges[,1] > 20000,
               ranges[, 2] - red.range[2] > 20000)

    disjoin.gr <- disjoin(subject.cnvs[j])
    disjoin.counts <- countOverlaps(disjoin.gr, subject.cnvs[j])

    par(mar = c(0, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
    layout(matrix(c(1,2,3,4,5,5), 2, 3, byrow=TRUE))

    ## Plot segments of chromosomal region first
    plot(NULL, xlim=c(min(ranges[,1]), max(ranges[,2])),
         ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position", ylab="Subject",
         xaxt="n", las=1, bty="l")
    axis(1, at=seq(c(min(ranges[,1]), max(ranges[,2]))), outer=TRUE)

    ## Plot line segments, color black if not in reduced range by over 20kb
    count <- 1
    for(l in 1:nrow(ranges)) {
        if(!w[l,1] && !w[l,2]){
            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2], y1=count,
                     col="gray")
        }
        else {
            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2],
                     y1=count, col="black")
        }
        count <- count + 1
    }
    abline(v=red.range, col="red", lty=2, lwd=2)

    ## Plot histogrom and fitted posterior
    posts <- getPost(r=r, kmin=1, kmax=5, delta=0.15,
                     S=2000, plot=T, burnin=500)
    text(max(r, na.rm=TRUE) - 0.1, max(hist(r, freq=FALSE)$density) - 0.1,
         paste("# Markers =", markers))
    ## Plot histogram and fitted posterior for adjusted log r ratios.
    posts1 <- getPost(r=r.adj, kmin=1, kmax=5, delta=0.15, S=2000,
                      plot=T, burnin=500)
    text(max(r.adj, na.rm=TRUE)-0.1, max(hist(r.adj, freq=FALSE)$density)-0.1,
         paste("# Markers =", markers))

    ## Plot density of overlapping regions
    xcoords <- c(rbind(start(disjoin.gr), end(disjoin.gr)))
    ycoords <- c(rbind(disjoin.counts, disjoin.counts))
    plot(NULL, xlim=c(min(x), max(y)), ylim=c(1, max(disjoin.counts)),
         las=1, bty="l")
    polygon(x=c(min(xcoords), xcoords, max(xcoords)), y=c(0, ycoords, 0),
            col=rgb(0,0.5,1, alpha=0.3), border="gray")
    abline(v=red.range, col="red", lty=2, lwd=2)

    ## Plot boxplots
    med <- median(df$avgRs)
    a <- aov(avgRs ~ plates, data=df)
    fstat <- summary(a)[[1]][["F value"]][1]
    graphics::boxplot(avgRs ~ plates, data=df, outline=FALSE, col="gray40",
                      names=zx, border="gray", xaxt="n", bty="n", ylim=c(-2.5, 1))
    axis(1, at=1:length(unique(zx)), las=2, cex.axis=0.4, labels=zx, tck=-0.01)
    abline(h=med, col="red", lty=2)
    if(fstat > 5) color="red"
    else color="black"
    text(1.2, -2.3, paste("F=",round(fstat, digits=2)), col=color)

    ## Name plot
    chr <- seqlevels(subject.cnvs[j])[table(seqnames(subject.cnvs[j])) != 0]
    mtext(text=paste(chr, ": ",region), side=3, outer=TRUE)
}
