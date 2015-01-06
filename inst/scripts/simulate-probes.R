library(Hmisc)
library(gtools)
## set seed

### Simulate probes
## l = number probes in each CNV
simulate.probes <- function(samples=2000, cnvs=200, K=4, probes=10, qual="easy") {
    l = probes

    Z <- array(data=NA, dim=c(cnvs,samples,K))
    Z[,,1] <- 1L
    ## Simulate copy numbers from multinomial distribution according to HWE
    theta = 30
    alpha = rbeta(200, 13.3, 13.3)
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
    }


    sl.good <- switch(qual,
                      easy = 6.25,
                      medium = 4.00,
                      hard = 3.00)
    #    sl.good <- c(6.25, 4.00, 3.00)
    sl.bad <- switch(qual,
                     easy = 0.0625,
                     medium = .400,
                     hard = .300)
    n <- switch(qual,
                easy = 0.2,
                medium = 0.5,
                hard = 0.5)

    ## Good quality
    ## dim = c(individual, probe, CNV, number of components)
    A <- array(data=NA, dim=c(samples, l, cnvs, K))
    Pb  <- matrix(NA, cnvs, l)
    Pv <- matrix(NA, cnvs, l)
    if(qual == "easy") for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, 0.03)
    if(qual == "medium" | qual == "hard")
        for(i in 1:nrow(Pb)) Pb[i,] <- rnorm(l, 0, 0.75)
    #    else stop("qual must be one of easy, medium, or hard")

    for(i in 1:nrow(Pb)) Pv[i,]   <- rgamma(l, shape = 19.92985, scale = 0.06272)

    ## medium and bad quality


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
easy.sim <- simulate.probes(samples = 8000, qual="easy")

saveRDS(easy.sim, file="easy-simulations.rds")
hist(rowMeans(easy.sim[[1]][,,15,3]), col="gray", breaks=80)

