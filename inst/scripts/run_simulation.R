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
model.files <- file.path(outdir, paste0("simulation_easy", seq_len(10 * 4), ".rds"))
##calledK <- matrix(NA, 10, 4)
##nIncorrect <- matrix(NA, 10, 4)
thr <- sapply(1:4, function(x) log(factorial(x)))
##thr[1] <- 0.5
##thr <- thr*3
thr <- c(0.5, 3, 6, 8)
mp <- McmcParams(iter=500, burnin=300, nStarts=10)
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
    x <- computeMarginalLik(pc, nchains=3,
                            T=1000, T2=500,
                            burnin=200,
                            K=1:4)
    saveRDS(x, file=model.files[it])
    if(FALSE){
      m <- readRDS(model.files[it])
      m2 <- m
      iter(m2) <- c(300, 600, 1000, 1000)
      nStarts(m2) <- 1
      burnin(m2) <- c(100, 400, 500, 500)
      m2 <- marginal(m2)
      saveRDS(m2, file=model.files[it])
      m <- m2
    }
    models <- orderModels(x)
    if(length(models) == 0){
      nc <- NA
    } else {
      nc <- k(models[[1]])
    }
    calledK[k] <- nc
  }
  ##list(calledK=calledK, nIncorrect=nIncorrect)
  return(calledK)
}
##simulation_summary <- list(K=calledK,
##                           nIncorrect=nIncorrect)
##saveRDS(simulation_summary,
##file="~/Software/CNPData/data/simulation_summary.rds")
dt <- Sys.Date()
saveRDS(results, file=paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))
q('no')

##results <- readRDS("~/Software/CNPData/data/simulation_summary.rds")
results <- readRDS("~/Software/CNPData/data/simulation_summary_maxperm3.rds")
calledK <- do.call(rbind, lapply(results, "[[", 1))
write.table(calledK, file="~/calledK2.txt", sep="\t")
write.csv(calledK, "")
truth <- matrix(1:4, nrow=10, ncol=4, byrow=TRUE)
nIncorrectCalls <- sum(truth != calledK)
## thresholds
## maxperm
#nIncorrect <- do.call(rbind, lapply(results, "[[", 2))


thr <- sapply(1:4, function(x) log(factorial(x)))
thr[1] <- 0.5
thr <- thr*3
calledK <- matrix(NA, 10, 4)
outdir <- "~/Software/CNPData/data"
model.files <- file.path(outdir, paste0("simulation", seq_len(10 * 4), ".rds"))
for(i in 1:10){
  for(k in 1:4){
    it <- (i-1)*4 + k
    m <- readRDS(model.files[it])
    s <- summary(m)
    s <- s[ s[, "range"] < 5, ]
    mns <- s[ , "mean"]
    ix <- order(mns, decreasing=TRUE)
    sel <- c(0, cumsum(abs(diff(mns[ix])) > 1))
    kk <- s[ix, "k"]
    kk <- min(kk[sel == 0])
    calledK[i, k] <- kk
##    m2 <- m
##    iter(m2) <- c(0, 300, 300, 300)
##    burnin(m2) <- c(0, 400, 400, 400)
##    nStarts(m2) <- 1
##    m2 <- marginal(m2)
##    m3 <- modelOtherModes(m2[[2]])
##
##    mp <- mcmcParams(m2[[2]])
##    burnin(mp) <- 0
##    mcmcParams(m2[[2]]) <- mp
##    sim <- posteriorSimulation(m3[[2]])
##    plot.ts(thetac(sim), plot.type="single", col=1:2)
##    bf <- bayesFactor(summary(m), thr=thr)
##    calledK[i, k] <- as.integer(substr(names(bf), 2, 2))
  }
}
k <- 1
i <- 2

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
