library(oligoClasses)
library(aricExperimentData)
library(aricUricAcid)
library(GenomicRanges)
library(snow)
library(doSNOW)
library(foreach)
NW <- 20
setParallelization(WORKERS=NW)

expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"
dat <- read.csv("/home/student/scristia/Software/CNPBayes/mccarroll-cnp.csv", header=TRUE, as.is=TRUE)

invisible(sapply(3:6, function(x) dat[,x] <<- as.integer(gsub(",", "", dat[,x]))))

ranges <- IRanges(dat[,5], dat[,6])
gr <- GRanges(seqnames=Rle(dat[,2]), ranges=ranges)
gr@elementMetadata@listData$CN <- dat[,7]

reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]
lrr.stats <- readRDS(file.path(uadir, "lrr_stats_penn.rds"))

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(chromosome(feature.gr) %in% autosomes())
feature.gr <- feature.gr[auto.index]

### chr1 regions we don't catch
#gr.chr1 <- gr[seqnames(gr)=="chr1"]
#reduced.chr1 <- reduced.gr[seqnames(reduced.gr)=="chr1"]
#hits <- subjectHits(findOverlaps(reduced.chr1, gr.chr1))

wc.files <- list.files(expdd, full.names=TRUE, pattern="wc_")
names(wc.files) <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
wc.files <- wc.files[as.character(pd.df$cel_id)] #8411

#set.seed(12)
#NR <- 15
#gr2 <- sample(gr.chr1[-hits], NR)
NR <- length(gr)

olaps <- findOverlaps(gr, feature.gr)
subj.index <- subjectHits(olaps)
query.index <- queryHits(olaps)


batchjob <- rep(seq_len(NW), length.out=length(wc.files))
fileList <- split(wc.files, batchjob)

## Find files that we need QC stats
fileList <- fileList[1:NW]
for(i in seq_len(NW)) fileList[[i]] <- fileList[[i]]
avgLrr <- foreach(k = seq_along(fileList), .combine="cbind") %dopar% {
    files <- fileList[[k]]
    A <- matrix(NA, NR, length(files))
    colnames(A) <- gsub(".rds", "", gsub("wc_", "", basename(files)))
    for(m in seq_along(files)){
        dat <- readRDS(files[m])
        rlist <- split(dat[subj.index, "robust_lrr"], query.index)
        A[, m] <- sapply(rlist, median, na.rm=TRUE)
    }
    A
}

se <- SummarizedExperiment(assays=SimpleList(avglrr=avgLrr),
                           rowData=gr)
#saveRDS(se, file="summarized_experiment_median_lrr_hg18.rds")
saveRDS(se, file="se_mccarroll.rds")
q('no')
