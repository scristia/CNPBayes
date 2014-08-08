library(cnvCall)
library(snow)
library(doSNOW)
library(foreach)
library(iterators)
library(devtools)
load_all("/home/bst/student/rscharpf/Software/HopkinsHPC/")
NW <- setParallelization(nc=20, force=TRUE)
avgRs <- readRDS("../../avgRs_wc_ea-vi.rds")
#lrr.list0 <- readRDS("~/MixtureModels/data/LRRList_adj.rds")
lrr.list0 <- readRDS("vi_lrr.rds")
ea <- match(rownames(avgRs), colnames(lrr.list0[[1]]))

nonas <- function(mat) {
    if(any(is.na(mat))) newmat <- mat[, -which(is.na(colSums(mat)==TRUE)),drop=FALSE]
    else newmat <- mat
    return(newmat)
}

## keep only EA with QC samples and transpose to be compatible with cnvCall
lrr.list1 <- sapply(lrr.list0, ea=ea, function(x, ea) { return(t(x[,ea,drop=FALSE])) } )
lrr.list2 <- sapply(lrr.list1, nonas)

arguments = set_arguments(write_results=F, quick=F, nocls=5:1, cnvs=1, collections="coll_1")
arguments$hier_mchoice_iters = 1:100
arguments$probe_summary      = "mean" # options: "pca","mean", "median"
arguments$signal_model       = "t8" # options: "t8", "normal"

PCAs <- matrix(NA, nrow(avgRs), ncol(avgRs))
rownames(PCAs) <- rownames(avgRs)
colnames(PCAs) <- paste0("VI_", seq_len(ncol(avgRs)))

batch <- sort(rep(1:NW, length.out=ncol(avgRs)))
batchjob <- split(1:length(batch), batch)

set.seed(4321)
#for(i in seq_along(lrr.list2)) {
pdf('~/tmp/cnvcall_pca-vi%03d.pdf', width=11, height=8, onefile=FALSE)
mat <- foreach(index=batchjob, .packages='cnvCall', .combine=cbind) %dopar% {
    A <- matrix(NA, nrow(avgRs), length(index))
    colnames(A) <- paste0("VI_",index)
    for(i in seq_along(index)) {
        data_matrix <- lrr.list2[[i]]
        qcallrv <- simple_call(data_matrix=data_matrix, simpleArguments = arguments)
        A[,i] <- qcallrv$ung_dat$coll_1
    }
    A
}
dev.off()
rownames(mat) <- rownames(avgRs)

saveRDS(mat, "PCA-vi.rds")

