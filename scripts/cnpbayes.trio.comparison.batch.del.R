#! /usr/bin/env Rscript

library(devtools)
#library(Rcpp)
library(tibble)
library(magrittr)
setwd("/dcs01/chaklab/chaklab1/users/mchou")
load_all("CNPBayes_trios")
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2")
ab <- commandArgs(trailingOnly=TRUE)  %>%
  as.integer
params.all <- readRDS(file.path("./params_del", paste0("params_", ab, ".rds")))
seeds <- readRDS(file.path("./params_del", "params_seeds.rds"))
seed <- seeds[ab]
set.seed(seed)

##--------------------------------------------------
##
message("cnpbayes_trios")
##
##--------------------------------------------------

# note maplabel and hp defined manually here for now

maplabel <- c(0,1,2)
hp <- HyperparametersTrios(k = 3)
N <- params.all[1,4]
params <- params.all[,1:3]
N<-1000
# sim truth
  nbatch <- 1
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0.001)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0.001, mprob, maplabel)
  
  is_offspring <- truth$data$family_member=="o"
  cn_offspring <- truth$data$copy_number[is_offspring]
  
  if (length(table(cn_offspring))!=hp@k) {
    truth <- simulate_data_multi2(params, N=N,
                                  batches = rep(c(1:nbatch),
                                                length.out = 3*N),
                                  error=0.001, mprob, maplabel)
  }
  
  # module for turning off ParamUpdates
  #up <- rep(1L, 10)
  #names(up) <- c("theta", "sigma2",
  #              "pi", "mu",
  #               "tau2",
  #              "nu.0",
  #               "sigma2.0",
  #               "z.parents",
  #               "z.offspring",
  #               "pi.parents")
  #up["z.parents"] <- 0L
  #up["z.offspring"] <- 1L
  
  #initialise model and mp
  #mp <- McmcParams(iter=4000, burnin=2000, thin=3, param_updates=up)
  #mp <- McmcParams(iter=4000, burnin=2000, thin=1)
  mp <- McmcParams(iter=3000, burnin=2000, thin=1)
  #model <- TBM(triodata=truth$data,
  #             hp=hp,
  #             mp=mp,
  #             mprob=mprob,
  #             maplabel=maplabel)
  
  
  for(i in 1:500){
    set.seed(123456)
    mp <- McmcParams(iter=1, burnin=i)
    model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                         batches=truth$data$batches,
                         mp=mp, k_range=c(3, 3), max_burnin=8000)[[1]]
    setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2/results")
    saveRDS(model, file=paste0("throwaway_model.",ab,".rds"))
  }
  
  #mp2 <- McmcParams(iter=4000, burnin=2000, thin=1)
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
               batches=truth$data$batches,
               mp=mp, k_range=c(3, 3), max_burnin=8000)

  #mcmcParams(model) <- mp
  #model <- posteriorSimulation(model)
  snr <- snr.calc(model[[1]])
  
  # this is for the trio model
  results <- z2cn(model[[1]], maplabel)
  is_offspring <- model[[1]]@triodata$family_member=="o"
  cn.truth <- truth$data$copy_number
  cn.truth.parent <- cn.truth[!is_offspring]
  cn.truth.child <- cn.truth[is_offspring]
  results.parent <- results@z[!is_offspring]
  results.child <- results@z[is_offspring]
  #cn.truth.parent2 <- as.factor(cn.truth.parent)
  #cn.truth.child2 <- as.factor(cn.truth.child)
  #results.parent2 <- as.factor(results.parent)
  #results.child2 <- as.factor(results.child)
  #cn.compare.parent <- table(results.parent2, cn.truth.parent2)
  #cn.compare.child <- table(results.child2, cn.truth.child2)
  
  prop.true.overall <- sum(cn.truth == results@z) / length(cn.truth)
  prop.true.parent <- sum(cn.truth.parent == results.parent) / length(cn.truth.parent)
  prop.true.child <- sum(cn.truth.child == results.child) / length(cn.truth.child)
  
  # these lines manually modified based on maplabel(del/dup)
  #sens.spec.par <- modified.sens.spec.calc(cn.compare.parent,cn.type="DEL")
  #sens.par <- sens.spec.par$sensitivity
  #spec.par <- sens.spec.par$specificity
  #sens.spec.chd <- modified.sens.spec.calc(cn.compare.child,cn.type="DEL")
  #sens.chd <- sens.spec.chd$sensitivity
  #spec.chd <- sens.spec.chd$specificity
  
  # this is for the MB model comparison
  results.mb <- z2cn(mb2[[1]], maplabel)
  #cn.truth2 <- as.factor(cn.truth)
  #results.mb2 <- as.factor(results.mb@z)
  #cn.compare.overall.mb <- table(results.mb2, cn.truth2)
  
  prop.true.overall.mb <- sum(cn.truth == results.mb@z) / length(cn.truth)
  
  #sens.spec.mb <- modified.sens.spec.calc(cn.compare.overall.mb,cn.type="DEL")
  #sens.mb <- sens.spec.mb$sensitivity
  #spec.mb <- sens.spec.mb$specificity
   true.stats <- component_stats(truth$data)
   truth.child.pi <- table(cn.truth.child)/length(cn.truth.child)
   child.pi <- table(results.child)/length(results.child)
   
   # input all the info
   summaryResults <- list(params = params.all,
                          SimTruth = true.stats,
                          SimChildPi = truth.child.pi,
                          TrioCNcall = results@z,
                          MBCNcall = results.mb@z,
                          SNR = snr, Accuracy = prop.true.overall,
                          AccuracyParents = prop.true.parent,
                          AccuracyChild = prop.true.child,
                          AccuracyMB = prop.true.overall.mb,
                          ModelPi = model[[1]]@pi,
                          ModelChildPi = child.pi,
                          ModelTheta = model[[1]]@theta,
                          ModelSigma2 = model[[1]]@sigma2
                          )
   

#epi.table <- cbind(SNR, accur.overall, accur.par, accur.child, accur.overall.mb)
#epi.table <- cbind(SNR, accur.overall, accur.par, accur.child, sens.parent, spec.parent, sens.child, spec.child, accur.overall.mb, sens.overall.mb, spec.overall.mb)
#names(epi.table) <- c("SNR", "ParentAccuracy", "ChildAccuracy", "ParSens", "ParSpec", "ChdSens", "ChdSpec")
#truth.params.df <- data.frame(truth.params.matrix)
#params.df <- data.frame(params.matrix)
#epi.table.df <- data.frame(epi.table)
#colnames(epi.table.df) <- c("SNR", "OverallAccuracy", "ParentAccuracy", "ChildAccuracy", "MBOverallAccuracy")
#colnames(epi.table.df) <- c("SNR", "OverallAccuracy", "ParentAccuracy", "ChildAccuracy", "ParSens", "ParSpec", "ChdSens", "ChdSpec", "MBOverallAccuracy", "MBSens", "MBSpec")

#summary.results <- vector("list", 3)
#summary.results[[1]] <- truth.params.df
#summary.results[[2]] <- params.df
#summary.results[[3]] <- epi.table.df

## save to some output directory
pathout <- file.path("/dcs01/chaklab/chaklab1/users/mchou/batchrun2/results")
#pathout <- file.path("~/Desktop/Chakravarti_Lab/git")
results.out <- paste0("params_", ab, ".txt")
#saveRDS(summary.results, file=file.path(pathout, params.rds))
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2/results")
write.table(as.data.frame(summaryResults),file=file.path(pathout, results.out), quote=F,sep="\t",row.names=F)
