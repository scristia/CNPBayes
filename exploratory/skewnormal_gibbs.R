## Skew normal mixture model
library(CNPBayes)
library(sn)
library(mvtnorm)
library(msm)
library(MASS)
library(mvtnorm)
library(gtools)
## real data
#avgRs <- readRDS("~/Projects/MixtureModel/avgRs.rds")
avgRs <- readRDS("../avgRs_wc_ea-vi.rds")
## 18, 347, 57, /!\78/!\, 107, 35
xx <- avgRs[,18]
xx <- xx[!is.na(xx)]
# simulated data (comment out when using real data)
omega <- c(4, 1)
omega2 <- omega^2
alpha <- c(-3, 0)
mu <- c(0, 4)

xx <- c(rsn(5000, mu[1], omega[1], alpha[1]), rsn(8000, mu[2], omega[2], alpha[2]))
xx <- xx[sample.int(8000)]
#plot(density(xx), type="l")
K <- 2
n <- length(xx)

##transformations
delta <- alpha/sqrt(1+alpha^2)
Ey <- mu+omega2*delta*sqrt(2/3.1415)
psi <- omega*delta
sigma2 <- omega2*(1-delta^2)


skewnormal.gibbs <- function(xx, priors, K, S, thin=1) {

    ### PRIORS

    beta0.psi <- 0
    beta0.xi <- mean(xx)
    D0.psi <- 0.1
    D0.xi <- 0.1
    B0 <- diag(c(D0.xi, D0.psi))
    phi <- 0.5
    c0 <- 2.5
    C0 <- phi * var(xx)
    #c0 <- 0.01
    #C0 <- 0.01
    eta <- rep(5, K)

    ## INITS
    #mu <- rep(mean(xx), K)
    alpha0 <- rep(0, K) ## skewness parameter
    #alpha0 <- c(-3, 0)
    omega0 <- rep(mad(xx), K) ## scale parameter
    omega20 <- omega0^2
    tau <- c(1,1)

    pi <- rep(1/K, K) ## intitial mixing params

    ## transformations
    delta <- alpha0/sqrt(1+alpha0^2)
    psi <- omega0*delta
    sigma20 <- omega20*(1-delta^2)
    tau <- 1/sigma20 ## precision

    beta0 <- rbind(rep(0,K), psi)

    ## starting values for Z and S: use kmeans, initialize Z to 0
    pars <- kmeans(xx, centers=K, nstart=15)
    mu <- sort(pars$centers)
    S <- rep(NA, length(xx))
    for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
    nn <- pars$size[order(pars$centers)]
    Z <- rep(0, length(xx))

    ### create storage matrices, initialize parameters from data
    nsim <- 2500
    thin <- 1
    thin.ind <- 1
    MU <- OMEGA <- ALPHA <- PREC <- PI <- matrix(rep(0, nsim/thin * K), ncol=K)
    assig <- matrix(0, nrow(xx), K)
    rownames(assig) <- rownames(xx)

    #MCMC
    for ( s in 1:nsim) {
        #draw Z
        v <- 1/(1+tau*psi^2)

        ## Update mu, alpha, omega
        for(k in 1:K) {
            ### Draw Z frum truncated normal distribution
#            v <- 1/(1+tau[k]*psi[k]^2)
            m <- v[k]*tau[k]*psi[k]*(xx[S==k]-mu[k])
            z <- rtnorm(nn[k], m, sqrt(v[k]), lower=0)

            ### don't use cbind, use matrix()
            X <- matrix(c(rep(1, length(z)), z), nrow=length(z))

            # Draw beta
            ### make these matrix operations faster
            ## diag(B0) <- 1/diag(B0)
#            B <- solve(B0 + tau[k]*crossprod(X))
#            beta <- B%*%(B0 %*% beta0[,k] + tau[k]*crossprod(X, xx[S==k]))
            B <- solve(B0 + tau[k]*crossprod(X))
            beta <- B%*%(B0 %*% beta0[,k] + tau[k]*crossprod(X, xx[S==k]))
            #    beta1 <- (xx[S==1] %*% X1 + c((1/D0.xi * beta0.xi), (1/D0.psi * beta0.psi)))%*%B1
            #    beta2 <- (xx[S==2] %*% X2 + c((1/D0.xi * beta0.xi), (1/D0.psi * beta0.psi)))%*%B2

            ## Draw xi and psi from their multivariate normal distribution
            #    mvdraw1 <- c(rmvnorm(1, c(beta1[1], beta1[2]), 1/tau[1] * B1))
            #    mvdraw2 <- c(rmvnorm(1, c(beta2[1], beta2[2]), 1/tau[2] * B2))
            mvdraw <- c(rmvnorm(1, beta, B))
            mu[k] <- mvdraw[1]
            psi[k] <- mvdraw[2]

            # Draw tau from its Gamma distribution
            cc <- c0 + nn[k]/2
#            eps <- crossprod(xx[S==k] - beta[2] - z * beta[k])
#            eps <- crossprod(xx[S==k] - X%*%mvdraw)
#            eps <- crossprod(xx[S==k] - X%*%beta)
#            C1 <- C0 + 0.5*(eps + 1/D0.xi*(beta[1] - beta0.xi)^2 +
#                            1/D0.psi*(beta[2] - beta0.psi)^2)

#            tau[k] <- rgamma(1, cc, C1)
            tau[k] <- rgamma(1, cc, C0 + crossprod(xx[S==k] - X%*%mvdraw)/2)
        }

        ## transformations
        alpha <- psi*sqrt(tau)
        omega <- sqrt(1/tau + psi^2)

        #### sample latent class variables
        d <- matrix(NA, nrow = length(xx), ncol = K)
        for(i in 1:K) d[,i] <- pi[i]*dsn(xx, mu[i], omega[i], alpha[i])
        p <- d/apply(d, 1, sum)
        #
        u <- runif(length((xx)))
        tmp <- p[, 1]
        S[u < tmp] <- 1
        if(K > 1){
            for(i in 2:K){
                S[tmp < u & u < tmp + p[, i]] <- i
                tmp <- tmp + p[, i]
            }
        }
        ##
        ## update [pi|data]
        ##
        for(i in 1:K) nn[i] <- sum(S==i)
        pi <- rdirichlet(1, nn+eta)

        # transform and store
        if(s%%thin == 0) {
            MU[thin.ind, ] <- mu
            OMEGA[thin.ind, ] <- omega
            ALPHA[thin.ind, ] <- alpha
            PI[thin.ind, ] <- pi
            thin.ind <- thin.ind+1
        }
        if(s%%100==0) cat(s,"\n")
    }
}

burnin <- 1:1000
mus <- colMeans(MU[-burnin, ])
omegas <- colMeans(OMEGA[-burnin, ])
alphas <- colMeans(ALPHA[-burnin, ])
pis <- colMeans(PI[-burnin, ])

source("sn_loglik.R")
snmixture <- list("mu"=MU[-burnin,], "omega"=OMEGA[-burnin,], "alpha"=ALPHA[-burnin,], "P"=PI[-burnin,])
loglik <- loglik.snmix(xx, snmixture, K)
bic <- -2*loglik + (4*K-1)*log(length(xx))
#9006.405 - 3, 9134.143 - 4, 


lim <- range(xx, na.rm=TRUE)
par(las=1, mfrow=c(1,1), mar=c(4, 4, 4, 4))
hist(xx, breaks = 500, col='lightgray', border='gray', freq=FALSE, main="",
     xlim=lim)
y2 <- seq(min(xx), max(xx), len=5000)
#post.dens <- pis[1]*dsn(y2, mus[1], omegas[1], alphas[1] ) + pis[2]*dsn(y2, mus[2], omegas[2], alphas[2])
#lines(y2, post.dens,lwd=2)
#mx <- max(post.dens)
for(k in 1:K) lines(y2, pi[k]*dsn(y2, mus[k], omegas[k], alphas[k] ), col="gray40", lty=2, lwd=2)
lines(y2, rowSums(sapply(1:K, function(x) pi[x]*dsn(y2, mus[x], omegas[x], alphas[x]))), col="skyblue3", lwd=2)
