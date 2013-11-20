# Check getPost parameters before running!`
library(GenomicRanges)
library(MixtureModel)
library(plyr)

## Read in genomic ranges data and avg R ratios
avgRs <- readRDS("~/MixtureModels/data/avgRs.rds")
granges <- readRDS("~/MixtureModels/data/granges.rds")
reduced <- readRDS("~/MixtureModels/data/reduced.rds")
grl <- readRDS("~/MixtureModels/data/vi_grl.rds")
subject.cnvs <- unlist(grl)

olap <- findOverlaps(reduced, subject.cnvs)
indices <- split(subjectHits(olap), queryHits(olap))

## find matching regions between granges and reduced regions
inds <- queryHits(findOverlaps(granges, reduced))

## READ IN PLATES
plates <- readRDS("~/MixtureModels/data/plate.rds")
lrr.list0 <- readRDS("~/MixtureModels/data/LRRList_red.rds")
nainds <- readRDS("~/nainds.rds")
samples <- readRDS("~/MixtureModels/data/samples.RDS")

source("~/MixtureModels/findnas.R")
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
w <- plates.rle$values[plates.rle$lengths < 50]

## names of plates with under 20 samples
wind <- which(!is.na(match(plates0, w)))
rmplates <- function(mat) mat[ ,-wind, drop=FALSE]

lrr.list2 <- sapply(lrr.list1, rmplates)

## average log2 ratios omitting plates with low number of samples
avgRs0 <- t(sapply(lrr.list2, colMeans))


## calculate f-statistics for each region using ANOVA
pval <- c()
fstat <- c()
for(i in seq_along(avgRs0)) {
    a <- aov(avgRs0[,i] ~ plates0[-wind])
    pval[i] <- summary(a)[[1]][["Pr(>F)"]][1]
    fstat[i] <- summary(a)[[1]][["F value"]][1]
}

### Center avgRs by median
r.list <- with(d, split(avgRs, plates))
meds <- sapply(r.list, median, na.rm=TRUE)
r.list.centered <- mapply(function(r, med) r-med, r=r.list, med=meds)
r.centered <- unlist(r.list.centered)


