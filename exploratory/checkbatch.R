# Check getPost parameters before running!`
library(GenomicRanges)
library(MixtureModel)
library(plyr)

## Read in genomic ranges data and avg R ratios
avgRs <- readRDS("~scristia/MixtureModels/data/avgRs.rds")
granges <- readRDS("~scristia/MixtureModels/data/granges.rds")
reduced <- readRDS("~scristia/MixtureModels/data/reduced.rds")
grl <- readRDS("~scristia/MixtureModels/data/vi_grl.rds")
subject.cnvs <- unlist(grl)

olap <- findOverlaps(reduced, subject.cnvs)
indices <- split(subjectHits(olap), queryHits(olap))

## find matching regions between granges and reduced regions
inds <- queryHits(findOverlaps(granges, reduced))

## READ IN PLATES
plates <- readRDS("~scristia/MixtureModels/data/plate.rds")
lrr.list0 <- readRDS("~scristia/MixtureModels/data/LRRList_adj.rds")
nainds <- readRDS("~scristia/nainds.rds")
samples <- readRDS("~scristia/MixtureModels/data/samples.RDS")

source("~scristia/MixtureModels/findnas.R")
## find and remove markers with NAs
lrr.list0 <- sapply(lrr.list0, function(x) x[, -nainds, drop=FALSE])
nonas <- function(mat) {
    mn <- !is.na(apply(mat, 1, sum))
    if(any(mn)) newmat <- mat[mn, ,drop=FALSE ]
    else newmat <- mn
    return(newmat)
}

lrr.list1 <- sapply(lrr.list0, nonas)
rm(nonas)

## Need to remove plates with less than 50 subjects
plates0 <- plates[-nainds]
plates.rle <- rle(plates0)
w <- plates.rle$values[plates.rle$lengths < 67]

## names of plates with under 20 samples
wind <- which(!is.na(match(plates0, w)))
rmplates <- function(mat) mat[ ,-wind, drop=FALSE]

lrr.list2 <- sapply(lrr.list1, rmplates)

## average log2 ratios omitting plates with low number of samples
avgRs0 <- t(sapply(lrr.list2, colMeans))


## calculate f-statistics for each region using ANOVA
pval <- rep(NA, length(avgRs0[,1]))
fstat <-rep(NA, length(avgRs0[,1]))
for(i in seq_along(avgRs0[,1])) {
    a <- aov(avgRs0[i,]/100 ~ plates0[-wind])
    pval[i] <- summary(a)[[1]][["Pr(>F)"]][1]
    fstat[i] <- summary(a)[[1]][["F value"]][1]
}

### Center avgRs by median
r.list <- with(d, split(avgRs, plates))
meds <- sapply(r.list, median, na.rm=TRUE)
r.list.centered <- mapply(function(r, med) r-med, r=r.list, med=meds)
r.centered <- unlist(r.list.centered)

library(lattice)
## exploratory by batch histograms for region 1
v <- rle(plates0[-wind])$lengths
vv <- sum(v[1:20])
vv2 <- sum(v[21:40])+vv
d <- data.frame("avgRs"=avgRs0[302,1:vv]/100, "plates"=plates0[-wind][1:vv])
d2 <- data.frame("avgRs"=avgRs0[302,(vv+1):vv2]/100, "plates"=plates0[-wind][(vv+1):vv2])
pdf("~/MixtureModels/batch/batch302_%01d.pdf", onefile=FALSE)
histogram(~avgRs | plates, d, breaks=40, border="gray", col="gray",
          panel=function(x,...) { panel.histogram(x,...)
          panel.abline(v=0, col="gray40", lty=2) } , main="Region 302")
histogram(~avgRs | plates, d2, breaks=40, border="gray", col="gray",
          panel=function(x,...) { panel.histogram(x,...)
          panel.abline(v=0, col="gray40", lty=2) } , main="Region 302")
dev.off()

