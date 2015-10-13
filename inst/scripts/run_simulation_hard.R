library(CNPBayes)
library(HopkinsHPC)
NW <- numberCores(10)
if(FALSE){
  arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                    "sl.bad" = 0.0625, ## sep param for "bad" probes
                    "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                    "n" = 0.2, ## background noise
                    "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                    )
  set.seed(123)
  dat <- CNPBayes:::simulateProbeLevel(cnvs=100, K=4, probes=10,
                                       arguments=arguments,
                                       qual="hard")
  saveRDS(dat, file="~/Software/CNPData/data/simulated_data_hard.rds")
} else {
  dat <- readRDS("~/Software/CNPData/data/simulated_data_hard.rds")
}

## dimensions are samples x probes x cnp x components
x <- dat[[1]]
cn <- dat[[2]]
outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation", seq_len(10 * 4), ".rds"))
mp <- McmcParams(iter=500, burnin=500, nStarts=10)
##mp <- McmcParams(iter=10, burnin=5, nStarts=3)

results <- foreach(i = 1:10, .packages=c("stats", "CNPBayes")) %dopar%{
  calledK <- rep(NA, 4)
  nIncorrect <- rep(NA, 4)
  cat(".")
  for(k in 1:4){
    it <- (i-1)*4 + k
    xx <- x[, , i, k]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
    if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
    model <- MarginalModel(data=pc)
    se <- as(model, "SummarizedExperiment")
    m <- marginal(se, mcmc.params=mp, maxperm=2, K=1:4)
    m2 <- m
    iter(m2) <- c(200, 500, 500, 500)
    ##iter(m2) <- iter(m2)/100
    nStarts(m2) <- 1
    ## need burnin for switching components
    burnin(m2) <- c(50, 300, 500, 750)
    ##burnin(m2) <- burnin(m2)/10
    m2 <- marginal(m2)
    saveRDS(m2, file=model.files[it])
    calledK[k] <- CNPBayes:::best(m)
    true_cn <- cn[i, , k]
    zz <- map(m[[k]])
    nIncorrect[k] <- sum(zz != true_cn)
  }
  list(calledK=calledK, nIncorrect=nIncorrect)
}
saveRDS(results, file="~/Software/CNPData/data/simulation_summary_hard_maxperm3.rds")

results <- readRDS("~/Software/CNPData/data/simulation_summary_hard_maxperm3.rds")
calledK <- do.call(rbind, lapply(results, "[[", 1))
nIncorrect <- do.call(rbind, lapply(results, "[[", 2))
