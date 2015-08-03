inits <- function(r, K, model="Normal"){
    ## NULL out these two variables to avoid NOTE about
    ## no visible binding for global variable
    ## SC: Can you take a look at this?
    xx <- nsim <- NULL
    if(model == "Normal") {
        if(K == 1) return(list("mu0"=mean(r), "sigma20"=var(r), "nn"=length(r)))
        ## Check if homozygous deletion may be present
        ## Also check greater than 1 observation so variance can be found.
        hdel <- min(r) < -0.9 & (length(r[r < -0.75]) > 1)
        ## if no: kmeans
        if(!hdel | sum(r < -0.75) > 1000) {
            pars <- kmeans(r, centers=K, nstart=50)
            mu0 <- sort(pars$centers)
            s20 <- pars$withinss/(pars$size - 1)[order(pars$centers)]
            nn <- pars$size[order(pars$centers)]
            return(list("mu0"=mu0,"sigma20"=s20,"nn"=nn))
        }
        ## if yes: kmeans (K-1) for data greater than -0.75
        else{
            mu1 <- median(r[r < -0.75])
            s201 <- 0.25^2
            pars <- kmeans(r[r > -0.75], centers=K-1, nstart=20)
            nn2 <- pars$size[order(pars$centers)]
            mu0 <- c(mu1, sort(pars$centers))
            s20 <- c(s201, pars$withinss/(pars$size - 1)[order(pars$centers)])
            nn <- c(length(r) - sum(nn2), nn2)
            return(list("mu0"=mu0,"sigma20"=s20,"nn"=nn))
        }
    }
    if(model == "SN") {
        ## initial values when skew normal mixture
        ## need: mu, omega (scale), alpha (skewness param), eta (mixing p)
        ## also: must return initial vector of allocations (from kmeans)
        if(K == 1) {
            alpha0 <- 0 ## skewness parameter
            omega0 <- mad(r) ## scale parameter
            omega20 <- omega0^2
            mu <- mean(r)
            S <- rep(0, length(xx))
            eta0 <- 1
            mat <- .Call("skewnormal_mix", r, K=K, S=S, centers=mu, alpha=alpha0,
                         omega2=omega20, eta=eta0, nsim)
            return(list("mu"=mean(r), "sigma20"=var(r), "nn"=length(r)))
        }
        ## Check if homozygous deletion may be present
        ## Also check greater than 1 observation so variance can be found.
        hdel <- min(r) < -0.9 & (length(r[r < -0.75]) > 1)
        ## if no: kmeans
        if(!hdel | sum(r < -0.75) > 1000) {
            pars <- kmeans(r, centers=K, nstart=50)
            mu0 <- sort(pars$centers)
            s20 <- pars$withinss/(pars$size - 1)[order(pars$centers)]
            nn <- pars$size[order(pars$centers)]
            return(list("mu0"=mu0,"sigma20"=s20,"nn"=nn))
        }
        ## if yes: kmeans (K-1) for data greater than -0.75
        else{
            mu1 <- median(r[r < -0.75])
            s201 <- 0.25^2
            pars <- kmeans(r[r > -0.75], centers=K-1, nstart=20)
            nn2 <- pars$size[order(pars$centers)]
            mu0 <- c(mu1, sort(pars$centers))
            s20 <- c(s201, pars$withinss/(pars$size - 1)[order(pars$centers)])
            nn <- c(length(r) - sum(nn2), nn2)
            return(list("mu0"=mu0,"sigma20"=s20,"nn"=nn))
        }
    }
}
#### Find components with significant overlap and call them as 'mixture of
#### mixtures'.


mix.dens <- function(y, comp, p, mix.mean, mix.sd) {
    p[comp]*dnorm(y, mean=mix.mean[comp], sd=mix.sd[comp])
}


## hack of code found on stackexchange, which I can not find again
min.f1f2 <- function(x, MU1, MU2, SD1, SD2, PI1, PI2) {
    f1 <- rowSums(sapply(1:length(PI1), mix.dens, y=x,
                          p=PI1, mix.mean=MU1, mix.sd=SD1), na.rm=TRUE)
    f2 <- PI2*dnorm(x, mean=MU2, sd=SD2)
    pmin(f1, f2)
}


### I don't use this
mergecomponent <- function(x, MU, SD, PI) {
    d <- diff(MU)/sqrt(sum(SD^2))
    area <- integrate(min.f1f2, -Inf, Inf, MU1=MU[1], MU2=MU[2], SD1=SD[1], SD2=SD[2], PI1=PI[1], PI2=PI[2])$value
    v <- max(area/PI)
    return(c(d, v))
}

mixture <- function(MU, SD, PI) {
    K <- length(MU)
    mixtures <- vector("list", K)
    i1 <- 1
    flag <- TRUE
    for(i in 1:(K-1)) {
        area <- integrate(min.f1f2, -Inf, Inf, MU1=MU[i1:i], MU2=MU[(i+1)],
                          SD1=SD[i1:i], SD2=SD[(i+1)], PI1=PI[i1:i], PI2=PI[(i+1)],
                          subdivisions=500L)$value
        v <- max(area/c(sum(PI[i1:i]), PI[(i+1)]))

        ### if area is greater than 0.5, treat as mixture
        if(v >= 0.5) {
            mixtures[[i]] <- i1:(i+1)
            if(length(i1:(i+1)) >= 3) mixtures[[i-1]] <- NULL
            flag <- TRUE
        }

        else if(v < 0.5) {
            if(flag == FALSE & i < K-1) {
                mixtures[[i]] <- i
                   }
            if(i < K-1 & flag == TRUE) {
                if(i == 1) mixtures[i] <- i
                mixtures[[i+1]] <-  i+1
                flag = FALSE
            }
            i1 <- i+1
            if(i == K-1) {
                if(flag == FALSE | i == 1) mixtures[[i]] <-  i
                mixtures[[i+1]] <-  i1
            }
        }
    }
    mixtures <- mixtures[!sapply(mixtures, is.null)]
    return(mixtures)
}
constr.draw <-
    function(mean, var, a, b){
        d <- pnorm(a, mean,
                   sqrt(var)) + runif(1) * (pnorm(b, mean, sqrt(var))
                   - pnorm(a, mean, sqrt(var)))
        theta <- qnorm(d, mean, sqrt(var))
        return(theta)
    }
dnormmix <-
    function(r, mixture, K, burnin=1, log=FALSE) {
        if(K == 1){
            mix.mean <- mean(mixture[["means"]][-burnin, ])
            mix.prec <- mean(mixture[["precs"]][-burnin, ])
            p <- mean(mixture[["P"]][-burnin, ])
        }
        else{
            mix.mean <- apply(mixture[["means"]][-burnin, ], 2, mean)
            mix.prec <- apply(mixture[["precs"]][-burnin, ], 2, mean)
            p <- apply(mixture[["P"]][-burnin, ], 2, mean)
        }
        ## Find likelihood for each component
        lik.comp <- function(r, comp) {
            p[comp]*dnorm(r, mean=mix.mean[comp], sd=1/sqrt(mix.prec[comp]))
        }
        ## Create likelihood array with likelihoods of from each component
        liks <- sapply(1:K, lik.comp, r=r)

        d <- rowSums(liks, na.rm=TRUE)
        if(log) {
            d <- log(d)
        }
        return(d)
    }
dppgibbs <- function(r, ## data
                     H, ## max number of clusters in block DP
                     alpha=1, ## concentration parameter
                     mu0, ## prior mean of theta for all clusters
                     tau20=0.1, ## prior prec for theta, all clusters
                     a=0.1,
                     b=0.1,
                     S=500
                     ){
    if(missing(H)) stop("H missing")
    n <- length(r) ## number of subjects
    ##
    #########
    # Inits #
    #########
    if(missing(mu0)){
        pars <- inits(r, H)
        mu0 <- mu <- pars$mu0
        tau20 <- tau <- 1/pars$sigma20
        sigma2 <- 1/tau
        ns <- pars$nn
        pi <- ns/n
    }

    #    pi<-ns<-rep(0,H)       # Mixing weights and number of subjects per cluster
    v<-rep(1/H,H)               # Conditional weights -- pr(c_i=h|c_i not in l<h)
    v[H]<-1                     # Apply DP truncation to H classes
    #    mu<-rep(0,H)           # Cluster-specific means
    #    tau<-sigma2<-rep(1,H)
    p<-tmp2<-matrix(0,n,H) # p[i,h] = posterior prob that subject i belongs to cluster h

    #########
    # Store #
    #########
    V<-Mu<-Sigma<-N<-Pi<-matrix(0,S,H)
    C<-matrix(0,S,n)
    grid<-seq(min(r),max(r),length=500)
    Y<-array(0,dim=c(S,length(grid),H))

    #########
    # GIBBS #
    #########
    for (i in 1:S) {
        # Update c, cluster indicator
        cumv<-cumprod(1-v)
        pi[1]<-v[1]
        for (h in 2:H) pi[h]<-v[h]*cumv[h-1]
        for (h in 1:H) tmp2[,h]<-pi[h]*dnorm(r,mu[h],sqrt(sigma2[h]))
        p<-tmp2/apply(tmp2,1,sum)
#        C[i,]<-c<-rMultinom(p,1)
        c<-rMultinom(p,1)
        Pi[i,]<-pi
        for (h in 1:H) ns[h]<-length(c[c==h])  # Must allow zeros for empty clusters

        # Update v
        for (h in 1:(H-1)) v[h]<-rbeta(1,1+ns[h],alpha+sum(ns[(h+1):H]))
        V[i,]<-v

        # Update mu and sigma2 and Yhat (density estimate)
        for (h in 1:H) {
            var<-1/(tau20+tau[h]*ns[h])
            m<-var*(tau20*mu0+tau[h]*sum(r[c==h]))
            Mu[i,h]<-mu[h]<-rnorm(1,m,sqrt(var))
            tau[h]<-rgamma(1,a+ns[h]/2,b+sum((r[c==h]-mu[h])^2)/2)
            Sigma[i,h]<-sigma2[h]<-1/tau[h]
            Y[i,,h]<-pi[h]*dnorm(grid,mu[h],sqrt(sigma2[h]))
        }
        N[i,]<-ns               # Number per cluster
        if (i%%100==0) print(i)

    }
    list(P=Pi, means=Mu, precs=1/Sigma, N=N)
}
dppgibbs <-
function(r, ## data
                     H, ## max number of clusters in block DP
                     alpha=1, ## concentration parameter
                     mu0=0, ## prior mean of theta for all clusters
                     tau20=0.1, ## prior prec for theta, all clusters
                     a=0.1,
                     b=0.1,
                     S=100
                     ){
        if(missing(H)) stop("H missing")
        n <- length(r) ## number of subjects
        ##
        #########
        # Inits #
        #########
         pi<-ns<-rep(0,H)       # Mixing weights and number of subjects per cluster
         v<-rep(1/H,H)          # Conditional weights -- pr(c_i=h|c_i not in l<h)
         v[H]<-1                        # Apply DP truncation to H classes
         mu<-rep(0,H)           # Cluster-specific means
         tau<-sigma2<-rep(1,H)
         p<-tmp2<-matrix(0,n,H) # p[i,h] = posterior prob that subject i belongs to cluster h

        #########
        # Store #
        #########
         V<-Mu<-Sigma<-N<-Pi<-matrix(0,S,H)
         C<-matrix(0,S,n)
         grid<-seq(min(y),max(y),length=500)
         Y<-array(0,dim=c(S,length(grid),H))

        #########
        # GIBBS #
        #########
        for (i in 1:S) {
         # Update c, cluster indicator
          cumv<-cumprod(1-v)
          pi[1]<-v[1]
          for (h in 2:H) pi[h]<-v[h]*cumv[h-1]
          for (h in 1:H) tmp2[,h]<-pi[h]*dnorm(y,mu[h],sqrt(sigma2[h]))
          p<-tmp2/apply(tmp2,1,sum)
          C[i,]<-c<-rMultinom(p,1)
          Pi[i,]<-pi
          for (h in 1:H) ns[h]<-length(c[c==h])  # Must allow zeros for empty clusters

         # Update v
          for (h in 1:(H-1)) v[h]<-rbeta(1,1+ns[h],alpha+sum(ns[(h+1):H]))
          V[i,]<-v

         # Update mu and sigma2 and Yhat (density estimate)
         for (h in 1:H) {
           var<-1/(tau20+tau[h]*ns[h])
           m<-var*(tau20*mu0+tau[h]*sum(y[c==h]))
           Mu[i,h]<-mu[h]<-rnorm(1,m,sqrt(var))
           tau[h]<-rgamma(1,a+ns[h]/2,b+sum((y[c==h]-mu[h])^2)/2)
           Sigma[i,h]<-sigma2[h]<-1/tau[h]
           Y[i,,h]<-pi[h]*dnorm(grid,mu[h],sqrt(sigma2[h]))
         }
         N[i,]<-ns              # Number per cluster
         if (i%%100==0) print(i)

        }
        list(P=Pi, means=Mu, precs=1/Sigma, Z=C, N=N)
}
## Convenience function for model selection and plotting
getPost <-
    function(r, tau20 = 0.1, nu0 = 1, kappa0=1, kmin=2, kmax=5, delta=0.15,
             S=200, plot=F, burnin=100, main="", crit="bic", full=FALSE){
        r <- r[!is.na(r)] ## remove missing values

        n <- length(r)
        loglik <- rep(NA, kmax - kmin + 1)
        bic <- rep(NA, kmax - kmin + 1)
        icl <- rep(NA, kmax - kmin + 1)
        kpost <- list()
        for(k in kmin:kmax){
            #nn <- rep(length(r)/k,k
            #mus <- kmeans(r, centers=k, nstart=15)$centers
            #mus <- sort(mus)
            inits <- inits(r,k)
            mu0 <- inits$mu0
            sigma20 <- inits$sigma20
            nn <- inits$nn
            print(nn)
            alpha <- rep(1,k)
            print(mu0)

            post <- gibbs.mix(r=r, k=k, tau20=tau20,
                              S=S, nu0=nu0, kappa0=1, burnin=burnin, delta=delta)

            loglik[k-kmin+1] <- post$loglik
            bic[k-kmin+1] <- post$bic
            icl[k-kmin+1] <- post$icl

            #save posterior
            kpost[[k-kmin+1]] <- post
        }
        ## print table showing log-likelihoods and BIC
        results <- data.frame("K"=seq(kmin,kmax), "log-likelihood"=loglik,
                              "BIC"=bic, "ICL"=icl)
        print(results)
        ## Return posterior of model with largest BIC
        if(plot){
            plotPosts(r=r, posts=kpost, burnin=burnin+1, main=main, crit=crit, full=full)
        }
        #if(!plot) return(kpost)
        if(crit == "icl")
            bestpost <- kpost[[which.min(sapply(1:length(kpost), function(x) kpost[[x]]$icl))]]
        else if(crit == "bic")
            bestpost <- kpost[[which.min(sapply(1:length(kpost), function(x) kpost[[x]]$bic))]]
        return(bestpost)
    }
gibbs <-
    function(r, nn,
             alpha=rep(1,3), ## dirichlet prior
             mu0=c(-4, -0.5, 0), ## theta ~ normal(mu0, tau20)
             tau20=0.1,  ## variance of theta
             nu0=1, ## number of prior observations for precision
             kappa0=1, ## number of prior observations for mu0
             sigma20=0.1, ## prec ~ gamma(nu0/2, nu0/2*sigma20)
             K=3, ## assume 3 components (for now)
             S=100,
             delta=0.15,
             burnin = 1:100
             ){

        is.1ornull <- function(x){(length(x) == 1) | is.null(x)}
        if(missing(mu0)){
            ## use kmeans for initial values
            pars <- inits(r, K)
            mu0 <- pars$mu0
            sigma20 <- pars$sigma20
            nn <- pars$nn
            print(mu0)
        }
        if(sum(nn) != length(r)) stop("nn must sum to length(r)")

        z <- rep(NA, length(r))
        s2n <- tau2n <- mun <- numeric(K)
        p <- matrix(NA, sum(nn), K)
        nun <- nu0+nn

        ## for sampling from constrained full conditionals
        a0 <- min(r)
        b0 <- max(r)

        s2 <- var(r, na.rm=TRUE)
        precs <- means <- matrix(NA, S, K)
        ## simulate from prior
        rbar <- means[1, ] <- rnorm(K, mu0, tau20)
        Z <- matrix(NA, length(r), S-1)
        ## just use marginal variance as guess of variance -- very diffuse
        #       precs[1,] <- 1/rep(s2, K)
        precs[1,] <- 1/sigma20
        PI <- matrix(NA,S, K)
        theta <- c()
        for(s in 2:S){
            ##
            ## simulate pi from its multinomial posterior
            ##
            pi <- rdirichlet(1, alpha+nn)
            ##
            ## update mu_n and tau_n^2
            ##
            tau2n <- 1/(1/tau20 + nn*precs[s-1, ])
            nun <- nu0+nn
            s2n <- 1/nun * (nu0*sigma20 + (nn-1)*s2 + kappa0*nn/(kappa0+nn)*(rbar - mu0)^2)
            mun <- (1/tau20)/(1/tau20 + nn*precs[s-1, ])*mu0 + nn*precs[s-1, ]/(1/tau20 + nn*precs[s-1, ])*rbar

            ##
            ## simulate from full conditional for theta
            ## samples from the constrained distributions for each theta
            ##
            ## This is not working properly.
            #           browser()
            theta <- means[s-1,]
            fint <- which(is.finite(theta))
            endpoints <- c(a0, theta[fint], b0)
            for(k in 1:length(fint)){
                if(is.finite(mun[fint[k]])){
                    a <- endpoints[k] + delta
                    b <- endpoints[k+2] - delta
                    if(mun[fint[k]] < a) mun[fint[k]] <- a
                    if(mun[fint[k]] > b) mun[fint[k]] <- b
                    theta[fint[k]] <- constr.draw(mun[fint[k]], tau2n[fint[k]], a, b)
                }
            }

            ##
            ## simulate precision from full conditional
            ##
            ## homozygous deletions should have large variance, so
            ## only keep prec for theta < 1 if std dev > 0.1
            prec <- rgamma(K, nun/2, nun/2 * s2n)
            #           homdel <- theta[is.finite(theta)] < -1
            #           if(any(homdel)){
            #                   if(any(1/sqrt(prec[which(homdel)]) < 0.1) & is.finite(prec[which(homdel)])){
            #                           prec[which(homdel)] <- precs[s-1, which(homdel)]
            #                   }
            #           }
            means[s, ] <- theta
            precs[s, ] <- prec

            ## simulate latent variable
            d <- matrix(NA, nrow = length(r), ncol = K)
            for(i in 1:K) d[,i] <- pi[i]*dnorm(r, theta[i], sqrt(1/prec[i]))
            p <- d/apply(d, 1, sum)

            #z <- rMultinom(p,1) - 1 ## Requires Hmisc package

            u <- runif(length(r))
            tmp <- p[, 1]
            z[u < tmp] <- 0
            if(K > 1){
                for(i in 2:K){
                    z[tmp < u & u < tmp + p[, i]] <- i-1
                    tmp <- tmp + p[, i]
                }
            }
            ##
            ## update [pi|data]
            ##
            for(i in 1:K) nn[i] <- sum(z==(i-1))

            ## Make sampling robust to when 1 observation is in a component
            ## or when 0 observations are in a component. Set as NA and only
            ## update remainder of the components.
            ##
            if(all(nn > 1)){
                rbar <- sapply(split(r, z), mean, na.rm=TRUE)
                s2 <- sapply(split(r, z), var, na.rm=TRUE)
            }
            else{
                ww <- split(r,z) ## Split by component membership
                ww <- ww[!sapply(ww, is.1ornull)] ## Remove components with 0 observations
                s2[nn > 0] <- sapply(ww, var, na.rm=TRUE)
                if(any(nn==1)) {
                    s2[nn == 1] <- NA
                    if(all(nn > 0)) rbar <- sapply(split(r,z), mean, na.rm=TRUE)
                }
                if(any(nn==0)){
                    s2[nn == 0] <- NA
                    rbar[nn != 0] <- sapply(split(r,z)[which(nn != 0)], mean, na.rm=TRUE)
                    rbar[nn == 0] <- NA
                }
            }

            ## for identifiability
            #           if(any(diff(rbar[is.finite(rbar)]) < 0)){
            #                   i <- which(diff(rbar[is.finite(rbar)]) < 0)
            #                   rbar[is.finite(rbar)][i+1] <- rbar[is.finite(rbar)][i]+0.01
            #           }

            Z[, s-1] <- z
            PI[s, ] <- pi
        }
        post <- list(P=PI, means=means, precs=precs, Z=Z, n=nn)
        loglik <- loglik.normmix(r, post, K=k, burnin=burnin)
        bic <- -2*loglik + (3*K-1)*log(length(r))

        c(post, "loglik"=loglik, "bic"=bic, "K"=K)
        #       list(P=PI, means=means, precs=precs, Z=Z, n=nn)
    }
gibbs.mix <- function(r, S=1000, k, delta=0.15, mu0, tau20,
                      nu0, sigma20, kappa0, burnin=100, outliers.rm=FALSE) {
    ## Check: if initial value vectors not of length k, STOP

  if(outliers.rm) {
        quant <- quantile(r, c(0.001, 0.999))
        rd <- r[r > quant[1] & r < quant[2] ]
        rd.low <- r[r <= quant[1]]
        rd.hi <- r[r >= quant[2]]
    }
    else rd <- r

    ## Check if initial values supplied (later)
    inits <- inits(rd, k)
    mu0 <- inits$mu0
#    sigma20 <- inits$sigma20
    nn <- inits$nn
    alpha <- rep(1, k)

    sigma20 <- rep(1, k)
    ## Initialize matrices for theta, prec, pi
    if(sum(nn) != length(rd)) stop("nn must sum to length(r)")

    z <- rep(NA, length(rd))
    s2n <- tau2n <- mun <- numeric(k)
    p <- matrix(NA, sum(nn), k)

    ## for sampling from constrained full conditionals
    a0 <- min(rd)
    b0 <- max(rd)

    ## For first iteration of mcmc, s2 is equal to sigma20
#    s2 <- sigma20
    s2 <- rep(mad(r), k)
    precs <- means <- PI <- matrix(NA, S, k)

    ## simulate from prior
#    rbar <- means[1, ] <- rnorm(k, mu0, tau20)
    rbar <- means[1, ] <- mu0
    Z <- matrix(0, length(rd), k, dimnames=list(names(rd), NULL))

    precs[1, ] <- 1/s2
    post <- .Call("gibbs_mix", r, means, precs, PI, Z, nu0, mu0, kappa0, alpha,
                  tau20, sigma20, rbar, s2, nn, delta, burnin)
#    post <- .Call("gibbs_mix_hier", r, means, precs, PI, Z, nu0, mu0, kappa0, alpha,
#                  tau20, sigma20, rbar, s2, nn, delta, burnin)

    if(outliers.rm) {
        post$Z <- rbind(post$Z, matrix(c(S-burnin, rep(0, k-1)),
                                       nrow=length(rd.low), ncol=k, byrow=TRUE,
                                       dimnames=list(names(rd.low), NULL)))
        post$Z <- rbind(post$Z, matrix(c(rep(0, k-1), S-burnin),
                                       nrow=length(rd.hi), ncol=k, byrow=TRUE,
                                       dimnames=list(names(rd.hi), NULL)))
    }

    MU <- colMeans(post$means)
    SD <- 1/sqrt(colMeans(post$precs))
    PI <- colMeans(post$P)
    CN <- getcn(MU)
    if(k > 1 & all(is.finite(SD) & is.finite(MU))) mixtures <- mixture(MU, SD, PI)
    else mixtures <- 1

    post <- c(post, list("CN"=CN))
    loglik <- loglik.normmix(rd, post, K=k, burnin=burnin)
    bic <- -2*loglik + (3*k-1)*log(length(rd))
    icl <- icl(rd, post, K=k, bic)

    return(c(post, "loglik"=loglik, "bic"=bic, "icl"=icl, "K"=k, "mix"=mixtures))
}

getcn <- function(mu) {
    if(any(is.na(mu))) return(NA)
    if(length(mu) == 1L) return(2)
    d <- diff(c(mu[1], mu[2]))
    if (mu < -0.7 || d > 0.6) {
        cn <- seq(0, length(mu)-1)
        return(cn)
    }
    else if (mu >= -0.7 && mu < -0.2 && d < 0.6) {
        cn <- seq(1, length(mu))
        return(cn)
    }
    else return(seq(2, length(mu)+1))
}
icl  <- function(r, post, K, bic) {
  ## find z_i,k (missing value variables)
  ## is.max: return vector of indicator variables, max(vec) = 1. If more than
  ## one max, choose one at random.
  is.max <- function(vec) {
    indic <- as.integer(vec == max(vec))
    ##if(sum(indic) > 1)
    invisible(indic)
  }
  z.hat<- t(apply(post$Z, 1, is.max))

  ## find t_k(x_i|theta_hat). This is conditional probability that x_i in
  ## component k given the parameter estimates.

  p.hat <- colMeans(post$P)
  theta.hat <- colMeans(post$means)
  sigma.hat <- 1/sqrt(colMeans(post$precs))
  ## tk is n dimensional vector of conditional probabilities each sample
  ## x_i belongs to a particular component.
  lik.comp <- function(r, comp) {
    p.hat[comp]*dnorm(r, mean=theta.hat[comp], sd=sigma.hat[comp])
  }

  tk <- sapply(1:K, lik.comp, r=r)
  tk <- tk/rowSums(tk)

  if(K!= 1) {
    if(all(tk > 0 & is.finite(tk))) entropy <- sum(colSums(tk * log(tk)))
    else {
      ## tentative until something better
      v <- tk*log(tk)
      v[which(!is.finite(v))] <- 0
      entropy <- sum(colSums(v))
    }
  }
  else entropy <- 0
  invisible(bic - 2*entropy)
}


setGeneric("ICL", function(object) standardGeneric("ICL"))

setMethod("ICL", "MarginalModel", function(object){


})
loglik.normmix <-
  function(r, mixture, K, burnin=1) {
    loglike <- dnormmix(r, mixture, K, burnin=1, log=TRUE)
    return(sum(loglike))
  }



logLikData <- function(object){
  B <- batch(object)
  mn <- theta(object)
  ss <- sigma(object)
  x <- y(object)
  tabz <- table(B, z(object))
  P <- tabz/rowSums(tabz)
  ## Vectorize
  lk <- k(object)
  xx <- rep(x, lk)
  nb <- rep(batchElements(object), lk)
  means <- rep(as.numeric(mn), nb)
  sds <- rep(as.numeric(ss), nb)
  p <- rep(as.numeric(P), nb)
  lik <- p*dnorm(xx, means, sds)
  lik <- matrix(lik, length(x), lk)
  lik <- rowSums(lik)
  sum(log(lik))
}

logLikData2 <- function(object){
  b <- batch(object)
  mn <- theta(object)
  ss <- sigma(object)
  x <- y(object)
  tabz <- table(b, z(object))
  P <- tabz/rowSums(tabz)
  B <- nBatch(object)
  K <- k(object)
  lik <- matrix(NA, length(b), K)
  for(i in 1:B){
    index <- b == i
    for(k in 1:K){
      lik[index, k] <-  P[i, k] * dnorm(x[index], mn[i, k], ss[i, k])
    }
  }
  lik <- rowSums(lik)
  sum(log(lik))
  ## Vectorize
##   lk <- k(object)
##   xx <- rep(x, lk)
##   nb <- rep(batchElements(object), lk)
##   means <- rep(as.numeric(mn), nb)
##   sds <- rep(as.numeric(ss), nb)
##   p <- rep(as.numeric(P), nb)
##   lik <- p*dnorm(xx, means, sds)
##   lik <- matrix(lik, length(x), lk)
##   lik <- rowSums(lik)
##   sum(log(lik))
}

logLikPhi <- function(object){
  thetas <- theta(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  sigma2s <- sigma2(object)
  ##
  ## Vectorize
  ##
  nr <- nrow(thetas)
  nc <- ncol(thetas)
  mus <- rep(mus, each=nr)
  taus <- rep(sqrt(tau2s), each=nr)
  thetas <- as.numeric(thetas)
  p.theta <- dnorm(thetas, mus, taus)
  p.theta <- matrix(p.theta, nr, nc)
  rownames(p.theta) <- uniqueBatch(object)
  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  sum(log(p.theta)) + sum(log(p.sigma2))
}

setMethod("computePrior", "BatchModel", function(object){
  compute_logprior_batch(object)
##   hypp <- hyperParams(object)
##   K <- k(hypp)
##   tau2s <- tau2(object)
##   mus <- mu(object)
##   p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
##   p.sigma2.0 <- dgamma(sigma2.0(object), shape=a(hypp), rate=b(hypp))
##   p.nu.0 <- dgeom(as.integer(nu.0(object)), betas(hypp))
##   sum(log(p.mu)) + log(p.sigma2.0) + log(p.nu.0)
})



.loglikMarginal <- function(object){
  x <- y(object)
  nr <- length(x)
  pp <- p(object)
  K <- k(object)
  x <- rep(x, K)
  p <- rep(p(object), each=nr)
  thetas <- rep(theta(object), each=nr)
  sigmas <- rep(sigma(object), each=nr)
  lik <- matrix(p*dnorm(x, thetas, sigmas),
                nr, K)
  sum(log(rowSums(lik)))
}

.loglikPhiMarginal <- function(object){
  thetas <- theta(object)
  mus <- mu(object)
  tau2s <- tau2(object)
  sigma2s <- sigma2(object)
  p.theta <- dnorm(thetas, mus, sqrt(tau2s))
  p.sigma2 <- dgamma(1/sigma2s, shape=1/2*nu.0(object), rate=1/2*nu.0(object)*sigma2.0(object))
  sum(log(p.theta)) + sum(log(p.sigma2))
}

.computeLoglik <- function(object){
  ##ll.data <- .loglikMarginal(object)
  ll.data <- loglik(object)
  ll.data
}

setMethod("computePrior", "MarginalModel", function(object){
  ##  .compute_prior_marginal(object)
  compute_logprior(object)
})

.compute_prior_marginal <- function(object){
  hypp <- hyperParams(object)
  K <- k(hypp)
  mus <- mu(object)
  p.sigma2.0 <- dgamma(sigma2.0(object), shape=a(hypp), rate=b(hypp))
  p.nu.0 <- dgeom(nu.0(object), betas(hypp))
  p.mu <- dnorm(mus, mu.0(hypp), sqrt(tau2.0(hypp)))
  log(p.mu) + log(p.sigma2.0) + log(p.nu.0)
}
## Skew normal mixture model
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

    ## starting values for mu, S and Z: use kmeans, initialize Z to 0
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
#            z <- constr.draw(m, v[k], a=0, b=Inf)
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
sntest <- function(r, K, nsim, burnin=500) {

    ##mu <- rep(mean(xx), K)
    xx <- r
    alpha0 <- rep(0, K) ## skewness parameter
    ##alpha0 <- c(-3, 0)
    omega0 <- rep(mad(xx), K) ## scale parameter
    omega20 <- omega0^2

    pars <- kmeans(xx, centers=K, nstart=15)
    mu <- sort(pars$centers)
    S <- rep(NA, length(xx))
    for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
    S <- S-1L
    nn <- pars$size[order(pars$centers)]

    eta0 <- rep(1/K, K) ## intitial mixing params
    mat <- .Call("skewnormal_mix", r, K=K, S=S, centers=mu, alpha=alpha0,
                 omega2=omega20, eta=eta0, nsim)

    return(mat)
}
## Plot genomic intervals across samples
## Input: subject.cnvs: subject specific copy number regions from HMM
##        red.range   : Reduced range across samples (from GenomicRanges)
##        indices     : Which samples to plot
##        region      : For graph (what genomic region is this?)
## Output:
##      Graph of all sample specific CNP regions as called by the HMM across
##      a region and a graph of the density of these overlapping regions.
##
## Todo: Add parameter option for whether to plot second graph.
plotSamples <- function(indices, red.range, subject.cnvs, region) {
    j <- indices
    x <- start(subject.cnvs)[j]
    y <- end(subject.cnvs)[j]
    ranges <- cbind(x, y)

    w <- cbind(red.range[1] - ranges[,1] > 20000,
               ranges[, 2] - red.range[2] > 20000)

    disjoin.gr <- disjoin(subject.cnvs[j])
    disjoin.counts <- countOverlaps(disjoin.gr, subject.cnvs[j])

    par(mar = c(0, 4, 0, 1), oma = c(4, 0, 4, 0) + 0.1, cex=0.5)
    layout(matrix(c(1,2,3,4,5,5), 2, 3, byrow=TRUE))

    ## Plot segments of chromosomal region first
    plot(NULL, xlim=c(min(ranges[,1]), max(ranges[,2])),
         ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position", ylab="Subject",
         xaxt="n", las=1, bty="l")
    axis(1, at=seq(c(min(ranges[,1]), max(ranges[,2]))), outer=TRUE)

    ## Plot line segments, color black if not in reduced range by over 20kb
    count <- 1
    for(l in 1:nrow(ranges)) {
        if(!w[l,1] && !w[l,2]){
            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2], y1=count,
                     col="gray")
        }
        else {
            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2],
                     y1=count, col="black")
        }
        count <- count + 1
    }
    abline(v=red.range, col="red", lty=2, lwd=2)

    ## Plot density of overlapping regions
    xcoords <- c(rbind(start(disjoin.gr), end(disjoin.gr)))
    ycoords <- c(rbind(disjoin.counts, disjoin.counts))
    plot(NULL, xlim=c(min(x), max(y)), ylim=c(1, max(disjoin.counts)),
         las=1, bty="l")
    polygon(x=c(min(xcoords), xcoords, max(xcoords)), y=c(0, ycoords, 0),
            col=rgb(0,0.5,1, alpha=0.3), border="gray")
    abline(v=red.range, col="red", lty=2, lwd=2)

    ## Name plot
    chr <- seqlevels(subject.cnvs[j])[table(seqnames(subject.cnvs[j])) != 0]
    mtext(text=paste(chr, ": ",region), side=3, outer=TRUE)
}
segplots <- function(subject.cnvs, red.range, indices, flag=20000L,
                     olap=FALSE, flag.col="gray", markers=NULL) {
    j <- indices ## which samples with CNVs
    x <- start(subject.cnvs)[j]
    y <- end(subject.cnvs)[j]
    ranges <- cbind(x, y)/1e6

    w <- cbind(start(red.range) - ranges[,1] > flag,
               ranges[, 2] - end(red.range) > flag)

    disjoin.gr <- disjoin(subject.cnvs[j])
    disjoin.counts <- countOverlaps(disjoin.gr, subject.cnvs[j])

    xlim <- c(max(0, start(red.range) - 1e4), end(red.range) + 1e4)/1e6
    if(sum(w[,1])/nrow(w) > 0.2) xlim[1] <- min(ranges[,1])
    if(sum(w[,2])/nrow(w) > 0.2) xlim[2] <- max(ranges[,2])
    ## Plot segments of chromosomal region first
    ## if overlapping == false, show x axis
    if(olap == FALSE) {
        plot(NULL, xlim=xlim, ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position",
             ylab="Subject", las=1, bty="l")
#        axis(1, at=seq(xlim[1], xlim[2], by = 5000L), outer=TRUE)
        ## show rug if markers provided
        if(!is.null(markers)) {
            ## rug for CN and SNP probes
            cn.ids <- grep("CN_", markers$feature.id)
            snp.ids <- grep("SNP_", markers$feature.id)
            cn.pos <- start(markers[cn.ids])/1e6
            snp.pos <- start(markers[snp.ids])/1e6

            rug(cn.pos, ticksize=-0.02, lwd=0.5)
            if(length(snp.pos) > 0)
                rug(snp.pos, ticksize=-0.03, lwd=1, col="blue")
        }
    }
    else {
        plot(NULL, xlim=xlim, ylim=c(0.7, nrow(ranges) + 0.3), xlab="Position",
             ylab="Subject", las=1, bty="l", xaxt='n')
    }
    #    axis(1, at=seq(xlim[1], xlim[2], by = 5000L), outer=TRUE)

    ## Plot line segments, color black if not in reduced range by threshold
    count <- 1
    sample.index <- seq_len(nrow(ranges))
    segments(x0=ranges[,1], y0=sample.index, x1=ranges[,2], y1=sample.index, col="gray")

#    for(l in seq_along(ranges[,1])) {
#        if(!w[l,1] && !w[l,2]){
#            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2], y1=count,
#                     col="gray")
#        }
#        else {
#            segments(x0=ranges[l,1], y0=count, x1=ranges[l,2],
#                     y1=count, col=flag.col)
#        }
#        count <- count + 1
#    }
    abline(v=c(start(red.range), end(red.range))/1e6, col="blue", lty=2, lwd=2)
    if(!olap & !is.null(markers)) {
        m2 <- markers[queryHits(findOverlaps(markers, red.range))]
        num.markers <- length(m2)
        legend('topleft', legend = paste0("Markers = ", num.markers), bty="n")
    }
    if(olap == TRUE) {
        xcoords <- c(rbind(start(disjoin.gr), end(disjoin.gr)))
        ycoords <- c(rbind(disjoin.counts, disjoin.counts))
        n <- max(disjoin.counts)
        overlapping(xcoords, ycoords, red.range, n, xlim, markers=markers)
    }

}

##
## Stephen, markers is not defined!  I added markers as an argument
## and pass it from the segplots function.
##
## need to fix coords to be on megabase scale
overlapping <- function(xcoords, ycoords, red.range, n, xlim, markers=NULL) {
    plot(NULL, xlim=xlim, ylim=c(1, n), las=1, bty="l")
    polygon(x=c(min(xcoords), xcoords, max(xcoords)), y=c(0, ycoords, 0),
            col=rgb(0,0.5,1, alpha=0.3), border="gray")
    abline(v=c(start(red.range), end(red.range)), col="blue", lty=2, lwd=2)

    if(!is.null(markers)) {
        ## rug for CN and SNP probes
        cn.ids <- grep("CN_", markers$feature.id)
        snp.ids <- grep("SNP_", markers$feature.id)
        cn.pos <- start(markers[cn.ids])
        snp.pos <- start(markers[snp.ids])

        num.markers <- length(snp.ids) + length(cn.ids)
        legend('topleft', legend = paste0("Markers = ", num.markers, bty='n'))

        rug(cn.pos, ticksize=-0.02, lwd=0.5)
        if(length(snp.pos) > 0)
            rug(snp.pos, ticksize=-0.04, lwd=1.5, col="blue")
    }
}
plotPosts <- function(r, posts, burnin=1, main="", crit="bic", full=TRUE) {
### This plotting function needs major cleaning up
## if a list of posterior output, do model selection and plot the best
    if(is.list(posts[[1]])) {
        logliks <- sapply(posts, function(x) x$loglik)
        K <- sapply(posts, function(x) x$K)
        bics <- sapply(posts, function(x) x$bic)
        icls <- sapply(posts, function(x) x$icl)
        dens.comp <- function(y, comp, p, mix.mean, mix.prec) {
            p[comp]*dnorm(y, mean=mix.mean[comp], sd=1/sqrt(mix.prec[comp]))
        }
#        y <- seq(min(r), max(r), len=1000)
        y <- seq(-3, 2, len=1000)
        ## Sort by BIC/ICL
        best.bic <- bics[order(bics)[1:2]]
        best.icl <- icls[order(icls)[1:2]]

        Kbic <- K[order(bics)[1:2]]
        if(crit == "bic") K2 <- K[order(bics)[1:2]]
        else if(crit == "icl") K2 <- K[order(icls)[1:2]]
        else stop("Only supported criterion are BIC and ICL")

        subs <- K2 - min(K) + 1
        ## Get posterior estimates of two best fitting models
        if(K2[1] == 1){
            p.1 <- mean(posts[[subs[1]]][["P"]])
            mean.1<- mean(posts[[subs[1]]][["means"]])
            prec.1 <- mean(posts[[subs[1]]][["precs"]])
        }
        else {
            p.1 <- apply(posts[[subs[1]]][["P"]], 2, mean)
            mean.1<- apply(posts[[subs[1]]][["means"]], 2, mean)
            prec.1 <- apply(posts[[subs[1]]][["precs"]], 2, mean)
        }
        ## second
        if(K2[2] == 1){
            p.2 <- mean(posts[[subs[2]]][["P"]][-burnin, ])
            mean.2<- mean(posts[[subs[2]]][["means"]])
            prec.2 <- mean(posts[[subs[2]]][["precs"]])
        }
        else{
            p.2 <- apply(posts[[subs[2]]][["P"]], 2, mean)
            mean.2<- apply(posts[[subs[2]]][["means"]], 2, mean)
            prec.2 <- apply(posts[[subs[2]]][["precs"]], 2, mean)
        }
        ## Draw histogram of data
        hist(r, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main=main, xlim=c(-3, 2))
        ## Overlay densities in different colors
        if(full) {
            d.1 <- rowSums(sapply(1:length(p.1), dens.comp, y=y,
                                  p=p.1, mix.mean=mean.1, mix.prec=prec.1), na.rm=TRUE)
            d.2 <- rowSums(sapply(1:length(p.2), dens.comp, y=y,
                                  p=p.2, mix.mean=mean.2, mix.prec=prec.2), na.rm=TRUE)
            lines(y, d.1, col="darkgreen", lwd=2)
            lines(y, d.2, col=rgb(24,167,181, maxColorValue=255), lty=2, lwd=2)
        }

        ###  plot individual components instead

        else{
            for(k in 1:K2[1]) {
                dens <- p.1[k]*dnorm(y, mean.1[k], sqrt(1/prec.1[k]))
                lines(y, dens, col='gray40', lwd=2)
                text(mean.1[k]-0.1, max(dens),
                     labels=posts[[subs[1]]]$CN[k], col="skyblue3")
#            for(k in 1:K2[1]) lines(y, p.1[k]*dnorm(y, mean.1[k], sqrt(1/prec.1[k])),
#                                    col='gray40', lwd=2)
            }
        }
        text(posts[[subs[1]]]$start, y=0, labels="x", col="tomato")
        ## Write key with K and BIC
        ## this is getting out of control
        legend("topleft",
               legend=paste(paste("K = ", c(paste(K2, collapse= " & "), paste(Kbic, collapse=" & ")) ,
                                  c(" ICL = "," BIC = "), c(paste(round(best.icl), collapse=" & "),
                                                            paste(round(best.bic), collapse=" & ")))), bty="n")
#        if(crit == "bic") {
#            legend("topleft", legend=paste(paste("K = ", K2, " BIC = ", round(best.bic))), col=c("darkgreen", rgb(24,167,181,max=255)), lty=1:2, bty="n")
#        }
#        else legend("topleft", legend=paste(paste("K = ", K2, " ICL = ", round(best.icl[1]))), col="gray40", lty=1, bty="n")
    }

## if posterior from a single model
    else {
        K <- posts$K
        y <- seq(-3, 2, len=1000)
        if(K == 1){
            p.1 <- mean(posts[["P"]])
            mean.1<- mean(posts[["means"]])
            prec.1 <- mean(posts[["precs"]])
        }
        else {
            p.1 <- apply(posts[["P"]], 2, mean)
            mean.1<- apply(posts[["means"]], 2, mean)
            prec.1 <- apply(posts[["precs"]], 2, mean)
        }
        ## Draw histogram of data
        hist(r, breaks=200, col="lightgray", border="lightgray", freq=FALSE,
             main=main, xlim=c(-3, 2))
        text(posts$start, y=0, labels="x", col="tomato")

        for(k in 1:K) {
            dens <- p.1[k]*dnorm(y, mean.1[k], sqrt(1/prec.1[k]))
            lines(y, dens, col='gray40', lwd=2)
            text(mean.1[k]-0.1, max(dens), labels=posts$CN[k], col="skyblue3")
        }
    }
}
