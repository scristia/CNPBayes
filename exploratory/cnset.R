library(Biobase)
library(GenomicRanges)
library(oligoClasses)
library(crlmm)
library(doSNOW)
library(snow)
library(ArrayTV)
library(ff)
library(devtools)
## if doesn't work, put inside foreach body
load_all("/home/bst/student/rscharpf/Software/ARIC/aricExperimentData")

setParallelization(WORKERS=15)
##outdir <- "/local_data/r00/aric/aricExperimentData/data"
outdir <- "/local/recover/r00/aric/aricExperimentData/data"
tmpdir <- "/local_data/r00/rscharpf/temp"
celdir <- "/dcs01/arking/ARIC_static/ARIC_Data/GWAS_CEL/freeze3_cel/"
expdd <- "/local/recover/r00/aric/aricExperimentData/data"
dir.create(tmpdir)
##
## fffiles are placed in a temporary directory that will be deleted
##
cels <- list.files(celdir, pattern=".CEL$", recursive=TRUE, full.names=TRUE)
plate <- substr(basename(cels), 1, 5)
cels.by.plate <- split(cels, plate)
cnset.names <- file.path(tmpdir, paste("cnSet_", names(cels.by.plate), ".rds", sep=""))

if (!exists('feature.gr'))
    feature.gr <- readRDS(file.path(expdd, "feature_granges.rds"))
autosomes <- paste0("chr",1:22)
auto.index <- which(seqnames(feature.gr) %in% autosomes)
feature.gr <- feature.gr[auto.index]

table(feature.gr$is.snp)
feature.snps <- feature.gr[feature.gr$is.snp]
### Store intensities all in big matrix instead
#bigA <- matrix(NA, length(feature.snps), length(cels), dimnames = list(feature.snps$feature.id, basename(cels)) )
#bigB <- matrix(NA, length(feature.snps), length(cels), dimnames = list(feature.snps$feature.id, basename(cels)) )

fid <- feature.snps$feature.id
foreach(i=seq_along(cnset.names),
        .packages=c("aricExperimentData", "oligoClasses",
                    "crlmm", "ff", "Biobase")) %dopar% {
    filename <- cnset.names[i]
    message("Processing ", filename)
    cels <- cels.by.plate[[i]]
    ##
    ## - should first check whether wave correction files exist.  If so, no reason to create cnSet
    ## - delete cnSets when all files present
    ldPath(tmpdir)
    ocProbesets(75000)
    ocSamples(100)
    celids <- cleanNames(basename(cels))
    ## create directory for ff objects
    cnSet <- aricExperimentData:::normalizeAB(cels, genome="hg18")
    cnSet <- cnSet[isSnp(cnSet), ]
    bigA <- A(cnSet)[fid,]
    bigB <- B(cnSet)[fid,]

    filenamesA <- file.path(outdir, paste0(celids,"_A.rds"))
    filenamesB <- file.path(outdir, paste0(celids,"_B.rds"))
    dimnames(bigA) = NULL
    dimnames(bigB) = NULL
    for(j in 1:ncol(cnSet)) {
       saveRDS(bigA[, j], file=filenamesA[j])
       saveRDS(bigB[, j], file=filenamesB[j])
       cat(".")
    }
}
