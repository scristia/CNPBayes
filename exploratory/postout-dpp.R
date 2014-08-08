library(snow)
library(doSNOW)
library(CNPBayes)
library(foreach)
library(iterators)
library(devtools)
load_all("/home/bst/student/rscharpf/Software/HopkinsHPC/")
NW <- setParallelization(nc=20, force=TRUE)
## VI regions
avgRs <- readRDS("../../avgRs_wc_ea-vi.rds")
colnames(avgRs) <- paste0("CNP", seq_along(avgRs[1,]))
## McCarroll regions
#avgRs <- readRDS("mc_avgRs_EA.rds")

S=1500
burnin=500
delta=1e-5
batch <- sort(rep(1:NW, length.out=ncol(avgRs)))
batchjob <- split(1:length(batch), batch)

set.seed(4321)
pst <- foreach(index=batchjob, .packages='CNPBayes') %dopar% {
    post <- vector("list", length(index))
    names(post) <- colnames(avgRs[,index])
    for(j in seq_along(index)) {
        k <- index[j]
        post[[j]] <- dppgibbs(avgRs[,k], H=5, alpha=1, S=S)
    }
    post
}
pst <- unlist(pst, recursive=FALSE)

saveRDS(pst, file="postdpp-vi.rds")
q('no')
