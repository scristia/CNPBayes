#! /usr/bin/env Rscript
library(devtools)
library(Rcpp)
library(tibble)
library(magrittr)
#setwd("/dcs01/chaklab/chaklab1/users/mchou")
##load_all("CNPBayes_trios")
svpacks(); setwd("trios")
load_all("CNPBayes")
set.seed(123456)
ab <- 10
##setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun2")
params.all <- readRDS(file.path("./params_dup", paste0("params_", ab, ".rds")))
maplabel <- c(2,3,4)
hp <- HyperparametersTrios(k = 3)
params <- params.all[,1:3]
N <- params.all[1,4]
n <- 5
SNR <- vector(length=n)
accur.overall <- vector(length=n)
accur.par <- vector(length=n)
accur.child <- vector(length=n)
accur.overall.mb <- vector(length=n)
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
i <- 1
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
mp <- McmcParams(iter=4000, burnin=2000, thin=3)
for(i in 100){
  set.seed(123456)
  mp <- McmcParams(iter=0, burnin=i)
  model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                       batches=truth$data$batches,
                       mp=mp, k_range=c(3, 3), max_burnin=8000)[[1]]
  saveRDS(model, file="throwaway_model.rds")
}
model <- readRDS("throwaway_models.rds")
set.seed(123456)
mp <- McmcParams(iter=0, burnin=1)
posteriorSimulation(model, mp) ## should fail with one iteration
## inspect inputs.  table(z(model)); theta(model), sigma(model),...


