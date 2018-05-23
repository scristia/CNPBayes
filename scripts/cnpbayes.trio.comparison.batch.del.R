#! /usr/bin/env Rscript

library(devtools)
library(Rcpp)
library(tibble)
library(magrittr)
setwd("/dcs01/chaklab/chaklab1/users/mchou")
load_all("CNPBayes_trios")
set.seed(123456)

ab <- commandArgs(trailingOnly=TRUE)  %>%
  as.integer
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2")

params.all <- readRDS(file.path("./params_del", paste0("params_", ab, ".rds")))

##--------------------------------------------------
##
message("cnpbayes_trios")
##
##--------------------------------------------------
#set.seed(193513)
## change as appropriate

#path <- "/dcs01/chaklab/chaklab1/users/mchou/batchrun2/params_del"
#params.all <- readRDS(file.path(path, params.rds))

# note maplabel and hp defined manually here for now
maplabel <- c(0,1,2)
hp <- HyperparametersTrios(k = 3)
params <- params.all[,1:3]
N <- params.all[1,4]

n <- 20
SNR <- vector(length=n)
accur.overall <- vector(length=n)
accur.par <- vector(length=n)
accur.child <- vector(length=n)
#sens.parent <- vector(length=n)
#sens.child <- vector(length=n)
#spec.parent <- vector(length=n)
#spec.child <- vector(length=n)
accur.overall.mb <- vector(length=n)
#sens.overall.mb <- vector(length=n)
#spec.overall.mb <- vector(length=n)

params.matrix <- matrix(nrow = n, ncol = 12, 
                        dimnames = list( c() , c("p1", "p2", "p3", 
                                                 "p1.child", "p2.child", "p3.child",
                                                 "theta1", "theta2", "theta3", 
                                                 "var1", "var2", "var3")))
truth.params.matrix <- matrix(nrow = n, ncol = 12, 
                              dimnames = list( c() , c("p1", "p2", "p3", 
                                                       "p1.child", "p2.child", "p3.child",
                                                       "theta1", "theta2", "theta3", 
                                                       "var1", "var2", "var3")))

for (i in 1:n) {

  nbatch <- 1
  #N <- 300
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
  
  #initialise model and mp
  mp <- McmcParams(iter=4000, burnin=2000, thin=3)
  
  #model <- TBM(triodata=truth$data,
  #             hp=hp,
  #             mp=mp,
  #             mprob=mprob,
  #             maplabel=maplabel)
  
  model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                       batches=truth$data$batches,
                       mp=mp, k_range=c(3, 3), max_burnin=8000)
  
  if (length(table(model[[1]]@z))!=hp@k){
    model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                         batches=truth$data$batches,
                         mp=mp, k_range=c(3, 3), max_burnin=8000)
  }
  
  mp2 <- McmcParams(iter=4000, burnin=2000, thin=1)
  mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
               batches=truth$data$batches,
               mp=mp2, k_range=c(3, 3), max_burnin=8000)
  
  if (length(table(mb2[[1]]@z))!=hp@k){
    mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
                 batches=truth$data$batches,
                 mp=mp2, k_range=c(3, 3), max_burnin=8000)
  }
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
  
  # input all the info
  SNR[i] <- snr
  accur.overall[i] <- prop.true.overall
  accur.par[i] <- prop.true.parent
  accur.child[i] <- prop.true.child
  #sens.parent[i] <- sens.par
  #spec.parent[i] <- spec.par
  #sens.child[i] <- sens.chd
  #spec.child[i] <- spec.chd
  accur.overall.mb[i] <- prop.true.overall.mb
  #sens.overall.mb[i] <- sens.mb
  #spec.overall.mb[i] <- spec.mb
  
  truth.child.pi <- table(cn.truth.child)/length(cn.truth.child)
  
  true.stats <- component_stats(truth$data)
  truth.params.matrix[i, 1:3] <- true.stats$p
  truth.params.matrix[i, 4:6] <- truth.child.pi
  truth.params.matrix[i, 7:9] <- true.stats$mean
  truth.params.matrix[i, 10:12] <- (true.stats$sd)^2
  
  # this section written specifically for k=3
  if (length(table(results.child))!=hp@k) {
    child.pi <- vector(length=hp@k)
    child.pi.tmp <- vector(length=3)
    child.pi.tmp[1] <- sum(results.child==maplabel[1])
    child.pi.tmp[2] <- sum(results.child==maplabel[2])
    child.pi.tmp[3] <- sum(results.child==maplabel[3])
    child.pi <- child.pi.tmp / length(results.child)
  } else {
    child.pi <- table(results.child)/length(results.child)
  }
  
  params.matrix[i, 1:3] <- model[[1]]@pi
  params.matrix[i, 4:6] <- child.pi
  params.matrix[i, 7:9] <- model[[1]]@theta
  params.matrix[i, 10:12] <- model[[1]]@sigma2
}

# remainder of lines for both loops

# to create ROC curve, otherwise ignore
#spec.parent2 <- 1-spec.parent
#spec.child2 <- 1 - spec.child

epi.table <- cbind(SNR, accur.overall, accur.par, accur.child, accur.overall.mb)
#epi.table <- cbind(SNR, accur.overall, accur.par, accur.child, sens.parent, spec.parent, sens.child, spec.child, accur.overall.mb, sens.overall.mb, spec.overall.mb)
#names(epi.table) <- c("SNR", "ParentAccuracy", "ChildAccuracy", "ParSens", "ParSpec", "ChdSens", "ChdSpec")
truth.params.df <- data.frame(truth.params.matrix)
params.df <- data.frame(params.matrix)
epi.table.df <- data.frame(epi.table)
colnames(epi.table.df) <- c("SNR", "OverallAccuracy", "ParentAccuracy", "ChildAccuracy", "MBOverallAccuracy")
#colnames(epi.table.df) <- c("SNR", "OverallAccuracy", "ParentAccuracy", "ChildAccuracy", "ParSens", "ParSpec", "ChdSens", "ChdSpec", "MBOverallAccuracy", "MBSens", "MBSpec")

summary.results <- vector("list", 3)
summary.results[[1]] <- truth.params.df
summary.results[[2]] <- params.df
summary.results[[3]] <- epi.table.df

## save to some output directory
pathout <- file.path("/dcs01/chaklab/chaklab1/users/mchou/batchrun2/results")
#pathout <- file.path("~/Desktop/Chakravarti_Lab/git")
results.out <- paste0("params_", ab, ".txt")
#saveRDS(summary.results, file=file.path(pathout, params.rds))
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2/results")
write.table(as.data.frame(summary.results),file=file.path(pathout, results.out), quote=F,sep="\t",row.names=F)
