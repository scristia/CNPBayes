
library(CNPBayes)
library(HopkinsHPC)
NC <- numberCores(10)
dat <- readRDS("~/Software/CNPData/data/simulated_data.rds")
## dimensions are samples x probes x cnp x components
x <- dat[[1]]
outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation_easy", seq_len(10 * 4), ".rds"))
##set.seed(1100)
T <- 1000; T2=500; burnin=250; thin <- 5; maxperm <- 3
index <- dim(x)[[3]]
results <- foreach(i = index, .packages=c("stats", "CNPBayes")) %dopar%{
  set.seed(i)
  calledK <- rep(NA, 4)
  nIncorrect <- rep(NA, 4)
  cat(i, " ")
  for(k in 1:4){
    it <- (i-1)*4 + k
    xx <- x[, , i, k]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center=TRUE, scale.=FALSE)$x[, 1]
    if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
    ds <- downSampleEachBatch(pc, nt=250, batch=rep(1, length(pc)))
    fit <- computeMarginalLik(ds$y,
                              nchains=3,
                              T=T, T2=T2,
                              burnin=burnin,
                              thin=thin,
                              maxperm=3,
                              K=1:4)
    model <- orderModels(fit)[[1]]
    dm <- DensityModel(model, merge=TRUE)
    calledK[k] <- k(dm)
  }
  return(calledK)
}
dt <- Sys.Date()
saveRDS(results, file=paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))
