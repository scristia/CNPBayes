library(oligoClasses)
library(CNPBayes)
library(aricExperimentData)
library(aricUricAcid)
library(GenomicRanges)
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
#gr <- sort(gr)

## using pennCNV calls
#reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')
reduced.gr <- readRDS(file.path(uadir, "reduced_granges.rds"))
vi.grl <- readRDS(file.path(uadir, "CNV_grl_EA_VI.rds"))
#vi.grl <- readRDS("/home/student/scristia/Software/CNPBayes/penncnv_granges.rds")
grv <- unlist(vi.grl)
dj <- disjoin(grv)

grv2 <- subsetByOverlaps(grv, reduced.gr)
j <- findOverlaps(grv2, reduced.gr, select="first")
grl <- split(grv2, j)

## setdiff do recursively
reduced.gr2 <- foreach(g=grl) %do%{
    ## g=grl[[1]]
    dj <- disjoin(g)
    cnt <- countOverlaps(dj, g)
    tmp <- reduce(dj[cnt > max(cnt)/2])##, min.gapwidth=1e3)

}
reduced.grl2 <- GRangesList(reduced.gr2)
reduced.gr2 <- unlist(reduced.grl2)

## load list of log r matrices
lrr <- readRDS("mccarroll_lrr.rds")
#lrr <- lrr[order(as.integer(names(lrr)))]
chv <-which(unlist(sapply(1:length(lrr), function(x) is.matrix(lrr[[x]]))) == FALSE)
lrr[[chv]] = matrix(lrr[[chv]], 1, length(lrr[[chv]]))

## number of probes
probes <- unlist(sapply(1:length(lrr), function(x) nrow(lrr[[x]])))
gr@elementMetadata@listData$probes <- probes

#ind <- unlist(sapply(1:length(lrr), function(x) nrow(lrr[[x]]) > 9))
#ind <- unlist(sapply(1:length(lrr), function(x) nrow(lrr[[x]]))) > 9
ind <- gr$probes > 5

lrr2 <- lrr[ind]

gr2 <- gr[as.integer(names(lrr2))]

hits <- subjectHits(findOverlaps(reduced.gr2, gr2))
gr3 <- gr2[-hits]

### phenotypes
pd.df <- readRDS(file.path(phenodir, "phenotypes.rds"))
pd.df <- pd.df[c(EA(pd.df), AA(pd.df)), ]
pd.df <- pd.df[keep(pd.df), ]

## EA samples
pd.ea <- pd.df[EA(pd.df), ]
pd.ea <- pd.ea[keep(pd.ea), ]

## AA samples
pd.aa <- pd.df[AA(pd.df), ]
pd.aa <- pd.aa[keep(pd.aa), ]
rm(pd.df)

### EA with quality controls
files <- list.files(expdd, pattern="vi_")
names(files) <- gsub(".rds", "", gsub("vi_", "", files))
files <- files[rownames(pd.ea)]

vi.grl <- readOrWriteFile(file.path(uadir, "CNV_grl_EA_VI.rds"),
                          collectGRanges,
                          datadir=expdd,
                          files=files,
                          max.cnv=100,
                          force=FALSE) ##8955

lrr.stats <- readRDS(file.path(expdd, "lrr_stats_penn.rds"))
lrr.stats <- lrr.stats[names(vi.grl), ]
drop <- which(elementLengths(vi.grl) > 100 | lrr.stats[, "lag10ACF"] > 0.03 | lrr.stats[, "mad"] > 0.32)
vi.grl <- vi.grl[-drop]
pd.ea <- pd.ea[names(vi.grl), ]
rm(lrr.stats, drop)

ea.names <- rownames(pd.ea)

ea.inds <- match(ea.names, colnames(lrr2[[1]]))

avgRs <- sapply(1:length(lrr2), function(x)
                colMeans(lrr2[[x]][,ea.inds, drop=FALSE]/100, na.rm=TRUE))
itx <- as.integer(names(lrr2))
colnames(avgRs) <- itx


#avgRs2 <- avgRs[,-hits]
#colnames(avgRs2) <- names(lrr2[-hits])

## mixture model params
delta=0.15
S=1000
burnin=200

fp <- 0
tp <- rep(NA, ncol(avgRs))
pdf('~/tmp/mccarroll-cnps%03d.pdf', width=11, height=8, onefile=FALSE)
jj <- 1
for(i in itx) {
    if(jj %% 4 == 1) par(mfrow=c(2,2), mar = c(4, 4, 4, 1), oma = c(4, 2, 4, 2) + 0.1, cex=0.5)
    chr <- chromosome(gr[i])
    main <- paste0(chr,": ", start(gr[i]), "-", end(gr[i]))
    tmp <- hist(avgRs[,paste(i)], plot=FALSE)
    set.seed(1234)
    posts <- getPost(r=avgRs[,paste(i)], kmin=1, kmax=5, delta=delta, S=S,
                     burnin=burnin, plot=TRUE, crit="icl", main=main)
    text(max(tmp$breaks), max(tmp$density), lab=gr$CN[i])
    K <- which.min(sapply(1:5, function(x) posts[[x]]$icl))
    if(K == 1) fp <- fp+1
    else if(K > 1) tp[i] <- i
    jj <- jj+1
}
dev.off()

tp <- tp[!is.na(tp)]
print(fp)
print(tp)

q('no')
