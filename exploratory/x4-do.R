library(snow)
library(doSNOW)
library(foreach)
library(GenomicRanges)
library(oligoClasses)
library(aricExperimentData)
library(aricUricAcid)
NW <- 20
setParallelization(WORKERS=NW)
expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]
lrr.stats <- readRDS(file.path(uadir, "lrr_stats_penn.rds"))

reduced.gr <- readRDS(file.path(uadir, "reduced_granges.rds"))

###########################################################
###########################################################

## using pennCNV calls
# vi.grl <- readRDS(file.path(uadir, "CNV_grl_EA_VI.rds"))
vi.grl <- readRDS("/home/student/scristia/Software/CNPBayes/penncnv_granges.rds")
gr <- unlist(vi.grl)
dj <- disjoin(gr)

gr2 <- subsetByOverlaps(gr, reduced.gr)
j <- findOverlaps(gr2, reduced.gr, select="first")
grl <- split(gr2, j)
## setdiff
reduced.gr2 <- foreach(g=grl) %do%{
    ## g=grl[[1]]
    dj <- disjoin(g)
    cnt <- countOverlaps(dj, g)
    reduce(dj[cnt > max(cnt)/2])##, min.gapwidth=1e3)
}
reduced.grl2 <- GRangesList(reduced.gr2)
reduced.gr2 <- unlist(reduced.grl2)
length(reduced.gr2) ## 378
median(width(reduced.gr2))##~13kb
NR <- length(reduced.gr2)

###########################################################
###########################################################


wc.files <- list.files(expdd, full.names=TRUE, pattern="wc_")
if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(chromosome(feature.gr) %in% autosomes())
feature.gr <- feature.gr[auto.index]

names(wc.files) <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
wc.files <- wc.files[as.character(pd.df$cel_id)] #8411

olaps <- findOverlaps(reduced.gr2, feature.gr)
subj.index <- subjectHits(olaps)
query.index <- queryHits(olaps)

batchjob <- rep(seq_len(NW), length.out=length(wc.files))
fileList <- split(wc.files, batchjob)

## Find files that we need QC stats
### Want a list of log Rs and BAF for all samples
fileList <- fileList[1:NW]
for(i in seq_len(NW)) fileList[[i]] <- fileList[[i]]
avgLrr <- foreach(k = seq_along(fileList), .combine="append") %dopar% {
    files <- fileList[[k]]
    A <- vector(mode="list", length = length(files))
    names(A) <- gsub(".rds", "", gsub("wc_", "", basename(files)))
#    for(m in seq_along(files)){
    for(m in 1:3) {
        dat <- readRDS(files[m])
        rlist <- split(dat[subj.index, "robust_lrr"], query.index)
        baflist <- split(dat[subj.index, "BAF"], query.index)
        ## create a matrix for each sample in list. col1: lrr, col2: baf
        A[[m]] <- mapply(lrr=rlist, baf=baflist, function(lrr,baf) cbind(lrr,baf))
    }
    A
}
se <- SummarizedExperiment(assays=SimpleList(avglrr=avgLrr),
                           rowData=reduced.gr2)
#saveRDS(se, file="summarized_experiment_median_lrr_hg18.rds")
saveRDS(se, file="se_median_lrr_penncnv.rds")
q('no')
