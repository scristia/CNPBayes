library(oligoClasses)
library(aricExperimentData)
library(aricUricAcid)
library(GenomicRanges)
library(foreach)

expdd <- "/local/recover/r00/aric/aricExperimentData/data"
uadir <- "/home/bst/student/rscharpf/Software/ARIC/aricUricAcid/data"
phenodir <- "/home/bst/student/rscharpf/Software/ARIC/aricPhenotypes/data"
dat <- read.csv("/home/student/scristia/Software/CNPBayes/mccarroll-cnp.csv", header=TRUE, as.is=TRUE)

invisible(sapply(3:6, function(x) dat[,x] <<- as.integer(gsub(",", "", dat[,x]))))

ranges <- IRanges(dat[,5], dat[,6])
gr <- GRanges(seqnames=Rle(dat[,2]), ranges=ranges)
gr@elementMetadata@listData$CN <- dat[,7]
autosomes <- paste0("chr",1:22)
gr <- gr[which(seqnames(gr) %in% autosomes)]
gr <- sort(gr)
reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
auto.index <- which(seqnames(feature.gr) %in% autosomes)
feature.gr <- feature.gr[auto.index]

olaps <- findOverlaps(gr, feature.gr)
subj.index <- subjectHits(olaps)
query.index <- queryHits(olaps)

pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]
wc.files <- list.files(expdd, full.names=TRUE, pattern="wc_")
names(wc.files) <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
wc.files <- wc.files[as.character(pd.df$cel_id)] #8411

nsamples <- length(wc.files)
samplenames <- gsub(".rds", "", gsub("wc_", "", basename(wc.files)))
M <- matrix(NA, length(subj.index), nsamples,
            dimnames=list(feature.gr$feature.id[subj.index], samplenames))

for(j in seq_len(nsamples)) {
    d <- readRDS(wc.files[j])
    M[,j] <- d[subj.index, "robust_lrr"]
}

indexList <- split(query.index, subj.index)
itx <- unique(unlist(indexList, use.names=FALSE))

matrixList <- foreach(i = itx) %do% M[which(query.index==i), ,drop==FALSE]
names(matrixList) <- itx

saveRDS(matrixList, file="mccarroll_lrr.rds")
q('no')
