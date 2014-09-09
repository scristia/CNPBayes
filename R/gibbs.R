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
        #	precs[1,] <- 1/rep(s2, K)
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
            #		browser()
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
            #		homdel <- theta[is.finite(theta)] < -1
            #		if(any(homdel)){
            #			if(any(1/sqrt(prec[which(homdel)]) < 0.1) & is.finite(prec[which(homdel)])){
            #			       	prec[which(homdel)] <- precs[s-1, which(homdel)]
            #			}
            #		}	
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
            #		if(any(diff(rbar[is.finite(rbar)]) < 0)){
            #			i <- which(diff(rbar[is.finite(rbar)]) < 0)
            #			rbar[is.finite(rbar)][i+1] <- rbar[is.finite(rbar)][i]+0.01
            #		}

            Z[, s-1] <- z
            PI[s, ] <- pi
        }
        post <- list(P=PI, means=means, precs=precs, Z=Z, n=nn)
        loglik <- loglik.normmix(r, post, K=k, burnin=burnin)
        bic <- -2*loglik + (3*K-1)*log(length(r))

        c(post, "loglik"=loglik, "bic"=bic, "K"=K)
        #	list(P=PI, means=means, precs=precs, Z=Z, n=nn)
    }
