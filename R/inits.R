inits <- function(r, K, model="Normal"){
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
