library(Hmisc)
library(gtools)
## set seed

arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                  "sl.bad" = 0.0625, ## sep param for "bad" probes
                  "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                  "n" = 0.2, ## background noise
                  "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                  )

sim <- function(samples=2000, cnvs=200, K=4, probes=10, arguments, qual="easy") {
    l <- probes
    sl.good <- arguments$sl.good
    sl.bad <- arguments$sl.bad
    prbias <- arguments$prbias
    n <- arguments$n
    prvar <- arguments$prvar

    Z <- array(data=NA, dim=c(cnvs,samples,K))
    Z[,,1] <- 1L
    ## Simulate copy numbers from multinomial distribution according to HWE
    theta = 30
    alpha = rbeta(200, 13.3, 13.3)
    ## probability of starting with homozygous,hemizygous,diploid
    ##
    for(j in 2:K) {
        k = j
        i = 0:(k-1)
        lambda = matrix(NA, cnvs, k)
        for(s in 1:cnvs) {
            w = choose(k-1, i) * alpha[s]^i * (1 - alpha[s])^((k-1) - i)
            lambda[s,] <- rdirichlet(1, w * theta)
        }
        #        w = choose(k-1, i) * alpha^i * (1 - alpha)^((k-1) - i)
        #        lambda <- rdirichlet(cnvs, w * theta)

        Z[,,j] <- rMultinom(lambda, samples)
        ##randomly shift multinomial sample to begin with 0, 1, or 2 copy number.
        ##Z[,,j] + sample(c(-1,0,1), 1)
    }

    A <- array(data=NA, dim=c(samples, l, cnvs, K))
    Pb  <- matrix(NA, cnvs, l)
    Pv <- matrix(NA, cnvs, l)
    for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, prbias)
    for(i in 1:nrow(Pb)) Pv[i,]   <- rgamma(l, shape = prvar[1], scale = prvar[2])

    corrupted <- sample(samples, ceil(0.004*samples))
    for(k in 1:K) {
        for(j in 1:cnvs) {
            for(i in 1:samples) {
                p.mean <- c(rep(sl.good,3),rep(sl.bad,7)) * Z[j,i,k] + Pb[j,]
                p.sd <- n + Pv[j,]

                A[i,,j,k] <- rnorm(l,  p.mean, p.sd)
            }
            ## null measurements for medium and bad quality
            if(qual == "medium" || qual == "hard") {
                v <- range(A[,,j,k])
                for(ind in corrupted) A[ind,,j,k] <- runif(l, v[1], v[2])
            }
        }
    }
    return(list("measurements"=A, "assignments"=Z))
}
### Simulate probes
## l = number probes in each CNV
simulate.probes <- function(samples=2000, cnvs=200, K=4, probes=10, qual="easy") {
    l = probes

    Z <- array(data=NA, dim=c(cnvs,samples,K))
    Z[,,1] <- 1L
    ## Simulate copy numbers from multinomial distribution according to HWE
    theta = 30
    alpha = rbeta(200, 13.3, 13.3)
    ## probability of starting with homozygous,hemizygous,diploid
    prob <- switch(qual,
                      easy = c(1/3,1/3,1/3),
                      medium = c(1/3,1/3,1/3),
                      hard = c(1/3,1/3,1/3))
    for(j in 2:K) {
        k = j
        i = 0:(k-1)
        lambda = matrix(NA, cnvs, k)
        for(s in 1:cnvs) {
            w = choose(k-1, i) * alpha[s]^i * (1 - alpha[s])^((k-1) - i)
            lambda[s,] <- rdirichlet(1, w * theta)
        }
        #        w = choose(k-1, i) * alpha^i * (1 - alpha)^((k-1) - i)
        #        lambda <- rdirichlet(cnvs, w * theta)

        Z[,,j] <- rMultinom(lambda, samples)
        ##randomly shift multinomial sample to begin with 0, 1, or 2 copy number.
        Z[,,j] + sample(c(-1,0,1),1,prob=prob)
    }


    ## controls separation of components, for good and bad probes.
    sl.good <- switch(qual,
                      easy = 6.25,
                      medium = 4.00,
                      hard = 3.00)
    sl.bad <- switch(qual,
                     easy = 0.0625,
                     medium = .400,
                     hard = .300)
    ## controls level of background noise
    n <- switch(qual,
                easy = 0.2,
                medium = 0.5,
                hard = 0.5)

    ## dim = c(individual, probe, CNV, number of components)
    A <- array(data=NA, dim=c(samples, l, cnvs, K))
    Pb  <- matrix(NA, cnvs, l)
    Pv <- matrix(NA, cnvs, l)
    if(qual == "easy") for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, 0.03)
    if(qual == "medium" | qual == "hard")
        for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, 0.75)
    #    else stop("qual must be one of easy, medium, or hard")

    for(i in 1:nrow(Pb)) Pv[i,]   <- rgamma(l, shape = 19.92985, scale = 0.06272)

    ## medium and hard


    #lth probe in ith individual
    corrupted <- sample(samples, ceil(0.004*samples))
    for(k in 1:K) {
        for(j in 1:cnvs) {
            for(i in 1:samples) {
                p.mean <- c(rep(sl.good,3),rep(sl.bad,7)) * Z[j,i,k] + Pb[j,]
                p.sd <- n + Pv[j,]

                A[i,,j,k] <- rnorm(l,  p.mean, p.sd)
            }
            ## null measurements for medium and bad quality
            if(qual == "medium" || qual == "hard") {
                v <- range(A[,,j,k])
                for(ind in corrupted) A[ind,,j,k] <- runif(l, v[1], v[2])
            }
        }
    }
    return(list("measurements"=A, "assignments"=Z))
}

set.seed(4321)
easy.sim <- simulate.probes(samples = 2000, qual="easy")

saveRDS(easy.sim, file="easy-simulations.rds")
hist(rowMeans(easy.sim[[1]][,,15,3]), col="gray", breaks=80)
hist(rowMeans(x[[1]][,,15,2])/10, col="gray", breaks=80)

