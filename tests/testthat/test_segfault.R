if(FALSE){
# before running the unit test, prep the folder structure and batch input files
# for reproducibility please use the params_seeds.rds file already in the testthat folder
# params_seeds.rds is currently for deletion seeds
#need to prep params tables in the form of cnp_1.rds

setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun6/params_dup")

paramsList <- vector("list", 12)

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[1]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[2]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 2000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[3]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[4]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[5]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 2000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[6]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[7]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[8]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 2000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[9]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[10]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[11]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 2000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[12]] <- params

# set up paramsList duplicates based on number of n replicates
n <- 10
paramsList2 <- rep(paramsList, n)
L <- length(paramsList2)

for(i in 1:L){
  params <- paramsList2[[i]]
  file<-paste("params_",i,".rds", sep="")
  saveRDS(params,file)
}

# make sure sample is sufficiently sized
# for reproducibility do not run
# seeds <- sample(1:1e6, replace=FALSE, size=5000)
# saveRDS(seeds, "params_seeds.rds")

# to reproduce the segfault
# set up directory with folder results and folder params_del
# put the params_seeds.rds file into the params_del folder
# put the deletion.script.R script in the main directory
# either run with shell script to run only particular input files
# or if running in qrsh, set the ab variable to following seeds:

#############
###rinse and repeat for deletions
###

setwd("/dcs01/chaklab/chaklab1/users/mchou/batchrun6/params_del")

paramsList <- vector("list", 30)

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[1]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 200
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[2]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[3]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 400
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[4]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[5]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 600
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[6]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 700
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[7]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 800
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[8]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 900
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[9]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[10]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[11]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 200
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[12]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[13]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 400
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[14]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[15]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 600
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[16]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 700
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[17]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 800
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[18]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 900
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[19]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[20]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[21]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 200
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[22]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[23]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 400
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[24]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[25]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 600
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[26]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 700
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[27]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 800
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[28]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 900
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[29]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[30]] <- params

# set up paramsList duplicates based on number of n replicates
n <- 10
paramsList2 <- rep(paramsList, n)
L <- length(paramsList2)

for(i in 1:L){
  params <- paramsList2[[i]]
  file<-paste("params_",i,".rds", sep="")
  saveRDS(params,file)
}

# make sure sample is sufficiently sized
# for reproducibility do not run
# seeds <- sample(1:1e6, replace=FALSE, size=5000)
# saveRDS(seeds, "params_seeds.rds")

# to reproduce the segfault
# set up directory with folder results and folder params_del
# put the params_seeds.rds file into the params_del folder
# put the deletion.script.R script in the main directory
# either run with shell script to run only particular input files
# or if running in qrsh, set the ab variable to following seeds:


}
