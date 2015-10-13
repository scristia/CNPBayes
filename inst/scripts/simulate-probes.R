library(Hmisc)
library(gtools)

## set seed
set.seed(4321)
## 'easy' arguments
arguments <- list("sl.good" = 6.25, ## separation parameter for "good" probes
                  "sl.bad" = 0.0625, ## sep param for "bad" probes
                  "prbias" = 0.03, ## probe level bias ~ N(0,prbias)
                  "n" = 0.2, ## background noise
                  "prvar" = c(19.92985, 0.06272) ## probe variance gamma parameters (shape,scale)
                  )

easy.sim <- CNPBayes:::simulateProbeLevel(samples=2000, cnvs=200, 
                                          arguments=arguments, qual="easy")

easy.sim <- simulate.probes(samples=2000, cnvs=200, arguments=arguments, qual="easy")

saveRDS(easy.sim, file="easy-simulations.rds")

hist(rowMeans(easy.sim[[1]][,,15,3]), col="gray", breaks=80)
hist(rowMeans(easy.sim[[1]][,,15,2])/10, col="gray", breaks=80)


### this moved to simulate_data.R
simulate.probes <- function(samples=2000, cnvs=200, K=4, probes=10, arguments, qual="easy") {
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
