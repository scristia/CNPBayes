## Skew normal mixture model
library(sn)
library(mvtnorm)
library(msm)
library(MASS)
library(mvtnorm)
library(gtools)
## real data
avgRs <- readRDS("~/Projects/MixtureModel/avgRs.rds")
xx <- avgRs[1, ]
xx <- xx[!is.na(xx)]
## simulated data (comment out when using real data)
#xx <- append(rsn(5000, 4, 1, 2.5), rsn(3000, 0, 1, 0))
K <- 2
phi <- 0.5
n <- length(xx)

### PRIORS
beta0.psi <- 0
beta0.xi <- mean(xx)
D0.psi <- 0.1
D0.xi <- 0.1
B0 <- diag(c(D0.xi, D0.psi))
c0 <- 2.5
C0 <- phi * var(xx)
#c0 <- 0.01
#C0 <- 0.01
eta <- rep(5, K)

## INITS
#mu <- rep(mean(xx), K)
alpha <- rep(0, K) ## skewness parameter
omega <- rep(var(xx), K) ## scale parameter
omega2 <- omega^2

pi <- rep(1/K, K) ## intitial mixing params

## transformations
delta <- alpha/sqrt(1+alpha^2)
psi <- omega*delta
sigma2 <- omega2*(1-delta^2)
tau <- 1/sigma2 ## precision

beta <- rbind(c(0,0), psi)

## starting values for Z and S: use kmeans, initialize Z to 0
pars <- kmeans(xx, centers=K, nstart=10)
mu <- sort(pars$centers)
S <- rep(NA, length(xx))
for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
nn <- pars$size[order(pars$centers)]
Z <- rep(0, length(xx))

### create storage matrices, initialize parameters from data
nsim <- 1000
thin <- 1
MU <- OMEGA <- ALPHA <- PREC <- PI <- matrix(rep(0, nsim/thin *K), ncol=K)

#MCMC
for ( s in 1:nsim) {
    ## update mixture probabilities
    pi <- rdirichlet(1, nn+eta)

    #draw Z
    v <- 1/(1+tau*psi^2)
    m <- v*tau*psi*(xx-mu)
    z1 <- rtnorm(nn[1], m[S==1], sqrt(v[1]), lower=0)
    z2 <- rtnorm(nn[2], m[S==2], sqrt(v[2]), lower=0)
    X1 <- cbind(rep(1, length(z1)), z1)
    X2 <- cbind(rep(1, length(z2)), z2)

    # Draw beta
    ### make these matrix operations faster
    B1 <- solve(solve(B0) + tau[1]*crossprod(X1))
    B2 <- solve(solve(B0) + tau[2]*crossprod(X2))
    beta1 <- B1%*%( solve(B0)%*%beta[,1] + tau[1]*crossprod(X1,xx[S==1]))
    beta2 <- B2%*%( solve(B0)%*%beta[,2] + tau[2]*crossprod(X2,xx[S==2]))
#    beta1 <- (xx[S==1] %*% X1 + c((1/D0.xi * beta0.xi), (1/D0.psi * beta0.psi)))%*%B1
#    beta2 <- (xx[S==2] %*% X2 + c((1/D0.xi * beta0.xi), (1/D0.psi * beta0.psi)))%*%B2

    ## Draw xi and psi from their multivariate normal distribution
#    mvdraw1 <- c(rmvnorm(1, c(beta1[1], beta1[2]), 1/tau[1] * B1))
#    mvdraw2 <- c(rmvnorm(1, c(beta2[1], beta2[2]), 1/tau[2] * B2))
    mvdraw1 <- c(rmvnorm(1, beta1, B1))
    mvdraw2 <- c(rmvnorm(1, beta2, B2))
    mu <- c(mvdraw1[1], mvdraw2[1])
    psi <- c(mvdraw1[1], mvdraw2[2])

    # Draw tau from its Gamma distribution
    cc <- c0 + nn/2
    eps1 <- crossprod(xx[S=1] - beta1[2] - z1 * beta1[1])
    eps2 <- crossprod(xx[S=2] - beta2[2] - z2 * beta2[1])
    #C1 <- C0 + 0.5*(eps1 + 1/D0.xi*(beta1[2] - beta0.xi)^2 + 1/D0.psi*(beta1[1] - beta0.psi)^2)
    #C2 <- C0 + 0.5*(eps2 + 1/D0.xi*(beta2[2] - beta0.xi)^2 + 1/D0.psi*(beta2[1] - beta0.psi)^2)
    tau1 <- rgamma(1, cc[1], C0 + crossprod(xx[S==1] - X1%*%mvdraw1)/2)
    tau2 <- rgamma(1, cc[2], C0 + crossprod(xx[S==2] - X2%*%mvdraw2)/2)
    tau <- c(tau1, tau2)

    ## transformations
    alpha <- psi*sqrt(tau)
    omega <- sqrt(1/tau + psi^2)

    #### sample latent class variables
    d <- matrix(NA, nrow = length(xx), ncol = K)
    for(i in 1:K) d[,i] <- pi[i]*dsn(xx, mu[i], omega[i], alpha[i])
    p <- d/apply(d, 1, sum)
#
#    S <- rMultinom(p, 1)
    u <- runif(length((xx)))
    tmp <- p[, 1]
    S[u < tmp] <- 1
    S[u >= tmp] <- 2
#    if(K > 1){
#        for(i in 2:K){
#            S[tmp < u & u < tmp + p[, i]] <- i
#            tmp <- tmp + p[, i]
#        }
#    }
    nn[1] <- sum(S==1)
    nn[2] <- sum(S==2)

    # transform and store
    if(s%%thin == 0) {
        MU[s, ] <- mu
        OMEGA[s, ] <- omega
        ALPHA[s, ] <- alpha
        PI[s, ] <- pi
    }
    if(s%%100==0) cat(s,"\n")
}


burnin <- 1:100
mus <- colMeans(MU[-burnin, ])
omegas <- colMeans(OMEGA[-burnin, ])
alphas <- colMeans(ALPHA[-burnin, ])
pis <- colMeans(PI[-burnin, ])

library(mcmcplots)
## need to name matrix columns before running this
MUS <- MU[-burnin, ]
OMEGAS <- OMEGA[-burnin, ]
ALPHAS <- ALPHA[-burnin, ]
PIS <- PI[-burnin, ]
colnames(MUS) <- c("mu1", "mu2")
colnames(PIS) <- c("pi1", "pi2")
colnames(OMEGAS) <- c("omega1", "omega2")
colnames(ALPHAS) <- c("alpha1", "alpha2")

mcmcplot(cbind(MUS, PIS, OMEGAS, ALPHAS), parms=c("mu1", "mu2", "pi1", "pi2",
                                             "omega1", "omega2", "alpha1", "alpha2"), dir="/Users/scrist/Dropbox/Mixture Models/skewnormal")

pdf("/Users/scrist/Dropbox/Mixture Models/skewnormal/sn_hist.pdf")
lim <- range(xx, na.rm=TRUE)
par(las=1, mfrow=c(1,1), mar=c(4, 4, 4, 4))
hist(xx, breaks = 500, col='gray', border='gray', freq=FALSE, main="",
     xlim=lim)
y2 <- seq(-3, 13, len=5000)
post.dens <- pis[1]*dsn(y2, mus[1], omegas[1], alphas[1] ) + pis[2]*dsn(y2, mus[2], omegas[2], alphas[2])
lines(y2, post.dens,lwd=2)
mx <- max(post.dens)
lines(y2, dsn(y2, mus[1], omegas[1], alphas[1] )/mx, col="gray40", lty=2)
lines(y2, dsn(y2, mus[2], omegas[2], alphas[2] )/mx, col="gray40", lty=2)
dev.off()



