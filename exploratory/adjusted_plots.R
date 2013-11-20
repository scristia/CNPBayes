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

## find indices of whether snp marker or CN
## (Not currently using these for anything)
cn.ind <- grep("CN", rownames(lrr.list0[[3]]))
snp.ind <- grep("SNP", rownames(lrr.list0[[3]]))

## Need to remove plates with less than 20 subjects
plates0 <- plates[-nainds]
plates.rle <- rle(plates0)
w <- plates.rle$values[plates.rle$lengths < 45]

## names of plates with under 20 samples
wind <- which(!is.na(match(plates0, w)))
rmplates <- function(mat) mat[ ,-wind, drop=FALSE]

lrr.list2 <- sapply(lrr.list1, rmplates)

## average log2 ratios omitting plates with low number of samples
avgRs0 <- t(sapply(lrr.list2, colMeans))

## Need list of F-statistics for each region, center regions with large 
## F values.
## Make one plot using rangePlots.R
source("~/MixtureModels/rangePlots.R")
avgRs.red  <- t(readRDS("~/MixtureModels/data/avgRs_red.rds"))
avgRs.adj  <- t(readRDS("~/MixtureModels/data/avgRs_adj.rds"))
gr.red     <- readRDS("~/MixtureModels/data/reduced.rds")
gr.adj     <- readRDS("~/MixtureModels/data/adjusted_granges.rds")
probes.red <- gr.red$numberProbes
probes.adj <- gr.adj$numberProbes

#pdf("~/MixtureModels/AllPlots/plots_09-24-13/adjustedplots%03d.pdf",
#    width=11, height=8, onefile=FALSE)
pdf("./testplot.pdf", width=11, height=8, onefile=FALSE)
set.seed(1234)
for(i in seq_along(reduced)) {
    df <- data.frame(plates = factor(plates0[-wind]), avgRs = avgRs0[i,])
    df$avgRs <- df$avgRs/100
    df <- arrange(df, plates)
    zx <- paste0(unique(df$plates), "(", table(df$plates), ")")
    region = paste0("Region ",i)

    red.range <- c(start(gr.red[i]), end(gr.red[i]))
    adj.range <- c(start(gr.adj[i]), end(gr.adj[i]))

    markers <- gr.red$numberProbes[i]
    markers.adj <- gr.adj$numberProbes[i]

    plotRegion(r = avgRs.red[i,], r.adj = avgRs.adj[i,], indices[[i]],
               red.range, adj.range, subject.cnvs, df, region, markers,
               delta=0.3, markers.adj, S=2200)
}
dev.off()
