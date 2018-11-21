############################################
### Use the same structure for batch prep
#############################################

# prep manifest file loop
# directory manifest
# reduce regions

setwd("/dcs01/chaklab/chaklab1/users/mchou/batchtest2/manifest")
for(i in 1:12){
  file<-paste("params_",i, sep="")
  write.table(NULL,file)
}

#need to prep params tables in the form of cnp_1.rds
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchtest2/params_del")

paramsList <- vector("list", 12)

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[1]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[2]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.05, 0.05, 0.05)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[3]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[4]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[5]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.1, 0.1, 0.1)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[6]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[7]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[8]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[9]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[10]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[11]] <- params

p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-1.5, 0.5)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[12]] <- params

for(i in 1:12){
  params <- paramsList[[i]]
  file<-paste("params_",i,".rds", sep="")
  saveRDS(params,file)
}

# duplications
# run once in qrsh
setwd("/dcs01/chaklab/chaklab1/users/mchou/batchtest2/params_dup")

paramsList <- vector("list", 12)

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[1]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[2]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.15, 0.15, 0.15)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[3]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[4]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[5]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.2, 0.2, 0.2)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[6]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[7]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[8]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.25, 0.25, 0.25)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[9]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 300
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[10]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 500
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[11]] <- params

p <- c(0.81, 0.18, 0.01)
theta <- c(0.5, 2, 4)
sigma2 <- c(0.3, 0.3, 0.3)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
paramsList[[12]] <- params

for(i in 1:12){
  params <- paramsList[[i]]
  file<-paste("params_",i,".rds", sep="")
  saveRDS(params,file)
}

