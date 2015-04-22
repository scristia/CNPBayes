## easy
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
  dat <- simulateProbeLevel(cnvs=100, K=4, probes=10,
                            arguments=arguments,
                            qual="easy")
  saveRDS(dat, file="~/Software/CNPData/data/simulated_data.rds")
} else {
  dat <- readRDS("~/Software/CNPData/data/simulated_data.rds")
}

## dimensions are samples x probes x cnp x components
x <- dat[[1]]
cn <- dat[[2]]
outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation", seq_len(10 * 4), ".rds"))
##calledK <- matrix(NA, 10, 4)
##nIncorrect <- matrix(NA, 10, 4)
mp <- McmcParams(iter=500, burnin=500, nStarts=10)

results <- foreach(i = 1:10, .packages=c("stats", "CNPBayes")) %dopar%{
  calledK <- rep(NA, 4)
  nIncorrect <- rep(NA, 4)
  cat(".")
  for(K in 1:4){
    it <- (i-1)*4 + K
    ##if(K < 3) next()
    ##    if(!file.exists(model.files[it])){
    if(TRUE){
      xx <- x[, , i, K]
      mns <- rowMeans(xx)
      pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
      if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
      model <- MarginalModel(data=pc)
      se <- as(model, "SummarizedExperiment")
      m <- marginal(se, mcmc.params=mp, maxperm=3, K=1:4)
      m2 <- m
      iter(m2) <- c(200, 500, 500, 500)
      nStarts(m2) <- 1
      ## need burnin for switching components
      burnin(m2) <- c(50, 300, 500, 750)
      saveRDS(m2, file=model.files[it])
    } else {
      ## if range for k=4 is in the 10-20 range, probably worth
      ## repeating with more burnin/iterations
      m <- readRDS(model.files[it])
      m2 <- m
      iter(m2) <- c(200, 500, 500, 500)
      nStarts(m2) <- 1
      burnin(m2) <- c(50, 300, 500, 750)
      m2 <- marginal(m2)
      saveRDS(m2, file=model.files[it])
      m <- m2
    }
    calledK[K] <- best(m)
    true_cn <- cn[i, , K]
    zz <- map(m[[K]])
    nIncorrect[K] <- sum(zz != true_cn)
  }
  list(calledK=calledK, nIncorrect=nIncorrect)
}
##simulation_summary <- list(K=calledK,
##                           nIncorrect=nIncorrect)
##saveRDS(simulation_summary,
##file="~/Software/CNPData/data/simulation_summary.rds")
saveRDS(results, file="~/Software/CNPData/data/simulation_summary_maxperm3.rds")

q('no')
##results <- readRDS("~/Software/CNPData/data/simulation_summary.rds")
results <- readRDS("~/Software/CNPData/data/simulation_summary_maxperm3.rds")
calledK <- do.call(rbind, lapply(results, "[[", 1))
nIncorrect <- do.call(rbind, lapply(results, "[[", 2))

outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation", seq_len(10 * 4), ".rds"))
k <- 4
i <- 6
it <- (i-1)*4 + k
m <- readRDS(model.files[it])
if(FALSE){
  plot(m[[1]])
  plot(m[[2]])
  plot(m[[3]])
  plot(m[[4]])
}
summary(m)
## when running again, use current values to start
m2 <- m
iter(m2) <- c(200, 500, 500, 500)
nStarts(m2) <- 1
## need burnin for switching components
burnin(m2) <- c(50, 300, 500, 750)
m2 <- marginal(m2)
summary(m2)

plot(m2[[4]], use.current=TRUE)

best(m2)
plot.ts(pic(m2[[4]]), plot.type="single", col=1:4)

calledK2 <- matrix(NA, 10, 4)
for(i in 1:10){
  for(k in 1:4){
    it <- (i-1)*4 + k
    m <- readRDS(model.files[it])
    calledK2[i, k] <- best(m)
  }
}
