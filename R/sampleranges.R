## Plot genomic intervals across samples
## Input: subject.cnvs: subject specific copy number regions from HMM
##        red.range   : Reduced range across samples (from GenomicRanges)
##        indices     : Which samples to plot
##        region      : For graph (what genomic region is this?)
## Output:
##      Graph of all sample specific CNP regions as called by the HMM across
##      a region and a graph of the density of these overlapping regions.
##
## Todo: Add parameter option for whether to plot second graph.
plotSamples <- function(indices, red.range, subject.cnvs, region) {
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

    ## Plot density of overlapping regions
    xcoords <- c(rbind(start(disjoin.gr), end(disjoin.gr)))
    ycoords <- c(rbind(disjoin.counts, disjoin.counts))
    plot(NULL, xlim=c(min(x), max(y)), ylim=c(1, max(disjoin.counts)),
         las=1, bty="l")
    polygon(x=c(min(xcoords), xcoords, max(xcoords)), y=c(0, ycoords, 0),
            col=rgb(0,0.5,1, alpha=0.3), border="gray")
    abline(v=red.range, col="red", lty=2, lwd=2)

    ## Name plot
    chr <- seqlevels(subject.cnvs[j])[table(seqnames(subject.cnvs[j])) != 0]
    mtext(text=paste(chr, ": ",region), side=3, outer=TRUE)
}
