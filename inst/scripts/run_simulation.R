library(CNPBayes)
library(HopkinsHPC)
NC <- numberCores(10)
if(FALSE){
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  set.seed(123)
  dat <- simulateProbeLevel(cnvs=100, K=4, probes=10,
                            arguments=arguments,
                            qual="easy")
  saveRDS(dat, file="~/Software/CNPData/data/simulated_data.rds")
} else {
  dat <- readRDS("~/Software/CNPData/data/simulated_data.rds")
}

## dimensions are samples x probes x cnp x components
x <- dat[[1]]
outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation_easy", seq_len(10 * 4), ".rds"))
set.seed(1100)
##results <- foreach(i = 1:10, .packages=c("stats", "CNPBayes")) %dopar%{
for(i in 1:10){
##  for(i in 5:10){
  calledK <- rep(NA, 4)
  nIncorrect <- rep(NA, 4)
  cat(".")
  for(k in 1:4){
    it <- (i-1)*4 + k
    xx <- x[, , i, k]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
    if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
    fit <- computeMarginalLik(pc, nchains=3,
                              T=2000, T2=1000,
                              burnin=1000,
                              K=1:4)
    models <- orderModels(fit)
    if(length(models) == 0){
      nc <- NA
    } else {
      nc <- k(models)[1]
    }
    calledK[k] <- nc
  }
##  return(calledK)
}
dt <- Sys.Date()
saveRDS(results, file=paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))
q('no')

dt <- Sys.Date()
results <- readRDS(paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))

##results <- readRDS("~/Software/CNPData/data/simulation_summary.rds")
results <- readRDS("~/Software/CNPData/data/simulation_summary_maxperm3.rds")
calledK <- do.call(rbind, lapply(results, "[[", 1))
write.table(calledK, file="~/calledK2.txt", sep="\t")
write.csv(calledK, "")
truth <- matrix(1:4, nrow=10, ncol=4, byrow=TRUE)
nIncorrectCalls <- sum(truth != calledK)
