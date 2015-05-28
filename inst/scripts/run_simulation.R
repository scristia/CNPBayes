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
##set.seed(1100)
T <- 2000; T2=1000; burnin=1000
results <- foreach(i = 1:10, .packages=c("stats", "CNPBayes")) %dopar%{
##for(i in 1:10){
  set.seed(i)
  calledK <- rep(NA, 4)
  nIncorrect <- rep(NA, 4)
  cat(i, " ")
  for(k in 1:4){
    it <- (i-1)*4 + k
    xx <- x[, , i, k]
    mns <- rowMeans(xx)
    pc <- prcomp(xx, center=TRUE, scale.=TRUE)$x[, 1]
    if(cor(pc, mns) < cor(-pc, mns)) pc <- -pc
    fit <- computeMarginalLik(pc, nchains=3,
                              T=T, T2=T2,
                              burnin=burnin,
                              K=1:4)
    models <- orderModels(fit)
    if(length(models) == 0){
      nc <- NA
    } else {
      nc <- k(models)[1]
    }
    calledK[k] <- nc
  }
  return(calledK)
}
dt <- Sys.Date()
saveRDS(results, file=paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))
q('no')

if(FALSE){
  library(cnvCall)
  data_matrix <- matrix(rnorm(1e5), ncol=10)
  data_matrix <- rbind(data_matrix, matrix(rnorm(1e5, mean=3, sd=1.4), ncol=10))
  rownames(data_matrix) <- paste("ind" , 1:nrow(data_matrix), sep="")
  colnames(data_matrix) <- paste("prob", 1:ncol(data_matrix), sep="")
  ##data_matrix <- matrix(as.numeric(data_matrix), ncol=1)
  ##dimnames(data_matrix) <- list(paste0("id", seq_len(nrow(data_matrix))), "cnp1")
  arguments <- set_arguments(write_results=F, quick=F, nocls=5:1, cnvs=1, collections="coll_1")
  arguments$hier_mchoice_iters <- 1:100
  arguments$probe_summary <- "pca" # or you could use: "t8""mean"
  arguments$signal_model <- "normal" # or you could use: "t8"
  qcallrv <- simple_call(data_matrix=data_matrix, simpleArguments = arguments)
}



dt <- Sys.Date()
dt <- "2015-05-12"
results <- readRDS(paste0("~/Software/CNPData/data/simulation_easy_", dt, ".rds"))

##results <- readRDS("~/Software/CNPData/data/simulation_summary.rds")
results <- readRDS("~/Software/CNPData/data/simulation_summary_maxperm3.rds")
calledK <- do.call(rbind, lapply(results, "[[", 1))
write.table(calledK, file="~/calledK2.txt", sep="\t")
write.csv(calledK, "")
truth <- matrix(1:4, nrow=10, ncol=4, byrow=TRUE)
nIncorrectCalls <- sum(truth != calledK)
