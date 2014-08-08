library(GenomicRanges)
library(CNPBayes)
m.gr <- readRDS(file="mccarroll-granges.rds")
vi.gr <- readRDS(file="VI-granges.rds")
posts <- readRDS(file="posteriors_all-vi.rds")
#posts <- readRDS(file="posteriors_EA.RDS")
mhet <- readRDS(file="baf-mat.rds")


getcn <- function(mu) {
    if(length(mu) == 1L) return(2)
    d <- diff(c(mu[1], mu[2]))
    mu1 <- abs(mu[1])
    if (mu1 > 0.7 | d > 0.6) {
        cn <- seq(0, length(mu)-1)
        return(cn)
    }
    else if (mu1 <=0.7 & mu1 < 0.2 & d < 0.6) {
        cn <- seq(1, length(mu))
        return(cn)
    }
    else return(seq(2, length(mu)+1))
}

mean1 <- sapply(1:length(posts), function(x) posts[[x]]$MU[1])
meandiff <- sapply(1:length(posts), function(x) diff(c(posts[[x]]$MU[1], posts[[x]]$MU[2])))
var1 <- sapply(1:length(posts), function(x) posts[[x]]$SIGMA[1])

CN <- sapply(1:length(posts), function(x) getcn(posts[[x]]$MU))
CN <- sapply(1:length(CN), function(x) paste(as.character(CN[[x]]), collapse=","))

vi.gr@elementMetadata$CN <- CN
vi.gr@elementMetadata$probes <- as.integer(NA)

m.gr@elementMetadata$source <- "McCarroll"
vi.gr@elementMetadata$source <- "VI"

hits <- findOverlaps(vi.gr, m.gr)

disc <- c(m.gr[-subjectHits(hits)], vi.gr[-queryHits(hits)])
disc <- sort(disc)

setGeneric("filterBy", function(query, subject, ...) standardGeneric("filterBy"))
setMethod("filterBy", c("GRanges", "GRanges"),
          function(query, subject, type="any", ...){
              dropbins <- unique(queryHits(findOverlaps(query, subject, type=type, ...)))
              if(length(dropbins) > 0)
                  query <- query[-dropbins]
              return(query)
          })
vi.gr2 <- filterBy(vi.gr, disc)
m.gr2 <- filterBy(m.gr, disc)

hits2 <- findOverlaps(vi.gr2, m.gr2)

vi.gr3 <- vi.gr2[queryHits(hits2)]
m.gr3 <- m.gr2[subjectHits(hits2)]

temp <- reduce(c(vi.gr3, m.gr3))
vi.regions <- findOverlaps(vi.gr3, temp, select="first")
vi.gr3$regions <- vi.regions

m.regions <- findOverlaps(m.gr3, temp, select="first")
m.gr3$regions <- m.regions

is.dupe <- duplicated(paste0(seqnames(vi.gr3), start(vi.gr3)))
vi.gr3 <- vi.gr3[!is.dupe]

is.dupe <- duplicated(paste0(seqnames(m.gr3), start(m.gr3)))
m.gr3 <- m.gr3[!is.dupe]

consensus  <- sort(c(vi.gr3, m.gr3))

saveRDS(consensus, paste0("cons-gr_", Sys.Date(), ".rds"))
saveRDS(disc, paste0("disc-gr_", Sys.Date(), ".rds"))



############
## need to think how to get component modes.
mean1 <- sapply(1:length(posts), function(x) MU[1])
meandiff <- sapply(1:length(posts), function(x) diff(c(MU[1], MU[2])))
var1 <- sapply(1:length(posts), function(x) SIGMA[1])
v <- cbind(mean1, meandiff, sqrt(var1))
#km <- kmeans(x=v[,c(1,2)], centers=3, nstart=15)$cluster
plot(v[,2], v[,1], pch=20, col="gray40")
abline(h=0.7)
abline(h=0.2)
abline(v=0.6)


assignments <- sapply(1:length(posts), function(x) apply(posts[[x]]$Z2, 1, which.max))
mhet <- mhet[match(rownames(assignments), rownames(mhet)), ]
ks <- sapply(1:length(posts), function(x) ncol(posts[[x]]$Z2))
assignments <- assignments[, ks>1]
mhet <- mhet[, ks>1]
v <- v[ks>2,]
cols <- rep("gray80", ncol(mhet))
cols[!is.nan(mhet[1,])] <- "gray20"
pdf("het-scplot.pdf")
plot(v[,2], v[,1], pch=20, col=cols, xlab="Mean difference", ylab="mean")
dev.off()

vvv <- vector(mode="list", length=ncol(mhet))
for(i in 1:ncol(mhet)) {
    vv <- split(mhet[,i]>0, assignments[,i])
    vvv[[i]] <- sapply(vv, mean)
}

posts <- readRDS("posts_vi2.rds")
avgRs <- readRDS("../../avgRs_wc_ea-vi.rds")
pdf('~/tmp/het-hists-bic%03d.pdf', width=11, height=8, onefile=FALSE)
for(i in seq_along(posts)) {
    plotPosts(avgRs[,i], posts[[i]], main=paste("Region", i))
    legend('topleft', legend =paste("comp", seq(1,length(vvv[[i]])),"=", round(vvv[[i]],3), collapse="\n"), bty="n")
}
dev.off()
