segplots <- function(subject.cnvs, red.range, indices, flag=20000L,
                     olap=FALSE, flag.col="gray", markers=NULL) {
    j <- indices ## which samples with CNVs
    x <- start(subject.cnvs)[j]
    y <- end(subject.cnvs)[j]
    ranges <- cbind(x, y)/1e6

    w <- cbind(start(red.range) - ranges[,1] > flag,
               ranges[, 2] - end(red.range) > flag)

    disjoin.gr <- disjoin(subject.cnvs[j])
    disjoin.counts <- countOverlaps(disjoin.gr, subject.cnvs[j])

    xlim <- c(max(0, start(red.range) - 1e4), end(red.range) + 1e4)/1e6
    if(sum(w[,1])/nrow(w) > 0.2) xlim[1] <- min(ranges[,1])
    if(sum(w[,2])/nrow(w) > 0.2) xlim[2] <- max(ranges[,2])
    ## Plot segments of chromosomal region first
    ## if overlapping == false, show x axis
    if(olap == FALSE) {
        plot(NULL, xlim=xlim, ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position",
             ylab="Subject", las=1, bty="l")
#        axis(1, at=seq(xlim[1], xlim[2], by = 5000L), outer=TRUE)
        ## show rug if markers provided
        if(!is.null(markers)) {
            ## rug for CN and SNP probes
            cn.ids <- grep("CN_", markers$feature.id)
            snp.ids <- grep("SNP_", markers$feature.id)
            cn.pos <- start(markers[cn.ids])/1e6
            snp.pos <- start(markers[snp.ids])/1e6

            rug(cn.pos, ticksize=-0.02, lwd=0.5)
            if(length(snp.pos) > 0)
                rug(snp.pos, ticksize=-0.03, lwd=1, col="blue")
        }
    }
    else {
        plot(NULL, xlim=xlim, ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position",
             ylab="Subject", las=1, bty="l", xaxt='n')
    }
    #    axis(1, at=seq(xlim[1], xlim[2], by = 5000L), outer=TRUE)

    ## Plot line segments, color black if not in reduced range by threshold
    count <- 1
    sample.index <- seq_len(nrow(ranges))
    segments(x0=ranges[,1], y0=sample.index, x1=ranges[,2], y1=sample.index, col="gray")

#    for(l in seq_along(ranges[,1])) {
#        if(!w[l,1] && !w[l,2]){
#            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2], y1=count,
#                     col="gray")
#        }
#        else {
#            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2],
#                     y1=count, col=flag.col)
#        }
#        count <- count + 1
#    }
    abline(v=c(start(red.range), end(red.range))/1e6, col="blue", lty=2, lwd=2)
    if(!olap & !is.null(markers)) {
        m2 <- markers[queryHits(findOverlaps(markers, red.range))]
        num.markers <- length(m2)
        legend('topleft', legend = paste0("Markers = ", num.markers), bty="n")
    }
    if(olap == TRUE) {
        xcoords <- c(rbind(start(disjoin.gr), end(disjoin.gr)))
        ycoords <- c(rbind(disjoin.counts, disjoin.counts))
        n <- max(disjoin.counts)
        overlapping(xcoords, ycoords, red.range, n, xlim)
    }

}

## need to fix coords to be on megabase scale
overlapping <- function(xcoords, ycoords, red.range, n, xlim) {
    plot(NULL, xlim=xlim, ylim=c(1, n), las=1, bty="l")
    polygon(x=c(min(xcoords), xcoords, max(xcoords)), y=c(0, ycoords, 0),
            col=rgb(0,0.5,1, alpha=0.3), border="gray")
    abline(v=c(start(red.range), end(red.range)), col="blue", lty=2, lwd=2)

    if(!is.null(markers)) {
        ## rug for CN and SNP probes
        cn.ids <- grep("CN_", markers$feature.id)
        snp.ids <- grep("SNP_", markers$feature.id)
        cn.pos <- start(markers[cn.ids])
        snp.pos <- start(markers[snp.ids])

        num.markers <- length(snp.ids) + length(cn.ids)
        legend('topleft', legend = paste0("Markers = ", num.markers, bty='n'))

        rug(cn.pos, ticksize=-0.02, lwd=0.5)
        if(length(snp.pos) > 0)
            rug(snp.pos, ticksize=-0.04, lwd=1.5, col="blue")
    }
}
