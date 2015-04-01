library(devtools)
load_all("~/Software/CNPBayes")
set.seed(100)
truth <- simulateData(N=2500,
                      theta=c(-2, -0.4, 0),
                      sds=c(0.3, 0.15, 0.15),
                      p=c(0.05, 0.1, 0.8))
if(FALSE) plot(truth, use.current=TRUE)
mcmcp <- McmcParams(iter=5000, thin=5, burnin=500, nStarts=20, nStartIter=200)
models <- fitMixtureModels(y(truth), mcmcp, K=1:4)
saveRDS(models, file="~/mm.rds")
q('no')

mcmcp <- McmcParams(iter=1000, thin=1, burnin=0)
models <- readRDS("~/mm.rds")
m3 <- models[[3]]
m3 <- useModes(m3)

partialGibbs <- function(model){
  m <- useModes(model)
  ptheta <- posteriorTheta(m, mcmcp)
  modes(m) <- modes(model)
  psigma <- posteriorSigma2(m, mcmcp)
  modes(m) <- modes(model)
  pmix <- posteriorP(m, mcmcp)
  results <- list(theta=ptheta,
                  sigma2=psigma,
                  p=pmix)
}
m4 <- permuteParameters(models[[4]])
par(mfrow=c(1,2), las=1)
m4 <- posteriorSimulation(m4, mcmcp)
m4 <- posteriorSimulation(models[[4]], mcmcp)
par(mfrow=c(1,2), las=1)
plot.ts(pic(m4), plot.type="single", col=1:4)
plot.ts(pic(models[[4]]), plot.type="single", col=1:4)




post3 <- partialGibbs(models[[3]])
post4 <- partialGibbs(models[[4]])
logp <- function(x) log(mean(x))
sapply(post3, logp)
sapply(post4, logp)

sum(sapply(post3, logp))
sum(sapply(post4, logp))
