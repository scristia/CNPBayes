library(oligoClasses)
library(CNPBayes)
library(aricExperimentData)
library(aricUricAcid)
library(aricPhenotypes)
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
gr <- sort(gr)
names(gr) <- paste0("MC", seq_along(gr))
#saveRDS(gr, "mc-granges.rds") ## without probe counts

## using VI CNV calls
reduced.gr <- readRDS("/home/student/scristia/Software/wcreduced.gr")
#reduced.gr <- readRDS('~/Software/wcreduced_penn.gr')

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
ind <- sort(gr)$probes > 5

lrr2 <- lrr[ind]

gr2 <- gr[gr$probes > 5]

#hits <- subjectHits(findOverlaps(reduced.gr, gr2))
hits <- subjectHits(findOverlaps(reduced.gr, gr))
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
aa.names <- rownames(pd.aa)

## get avgRs for EA and AA samples
ea.inds <- match(ea.names, colnames(lrr2[[1]]))
aa.inds <- match(aa.names, colnames(lrr2[[1]]))


avgRs <- sapply(1:length(lrr), function(x)
                colMeans(lrr[[x]][,c(ea.inds, aa.inds),
                         drop=FALSE]/100, na.rm=TRUE))

## make sure region is from correct location in granges object
itx <- as.integer(names(lrr))
colnames(avgRs) <- itx


#avgRs2 <- avgRs[,-hits]
#colnames(avgRs2) <- names(lrr2[-hits])

## mixture model params
delta=0.15
S=1000
burnin=200

pdf('/home/student/scristia/tmp/mccarroll-cnps%03d.pdf', width=11, height=8, onefile=FALSE)
jj <- 1
posteriors <- list()
length(posteriors) <- ncol(avgRs)
for(i in itx) {
    if(jj %% 4 == 1) par(mfrow=c(2,2), mar = c(4, 4, 4, 1), oma = c(4, 2, 4, 2) + 0.1, cex=0.5)
    chr <- chromosome(gr[i])
    main <- paste0(chr,": ", start(gr[i]), "-", end(gr[i]))
    tmp <- hist(avgRs[,paste(i)], plot=FALSE)
    set.seed(1234)
    posts <- getPost(r=avgRs[,paste(i)], kmin=1, kmax=5, delta=delta, S=S,
                     burnin=burnin, plot=TRUE, crit="icl", main=main)
    text(max(tmp$breaks), max(tmp$density), lab=gr$CN[i])
    K <- posts$K
    PI <- colMeans(posts$P)
    MU <- colMeans(posts$means)
    SIGMA <- colMeans(1/posts$precs)
    Z <- posts$Z

    ## need to find markers
#    snp.ids <- grep("SNP_", m$feature.id)
#    snp.pos <- start(m[snp.ids])
#    posteriors[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, snps=snp.pos, gr=gr[i])
    posteriors[[i]] <- list(K=K, PI=PI, MU=MU, SIGMA=SIGMA, Z=Z, gr=gr[i])
    jj <- jj+1
}
dev.off()
saveRDS(posteriors, "mc-posteriors.rds")

q('no')


pr5 <- sapply(1:length(xlist2), function(x) xlist[[x]]$gr$probes) > 5
xlist3 <- xlist2[pr5]
length(xlist2[pr5])

saveRDS(xlist3, "mc-posteriors2.rds")
