dsnmix <-
    function(r, mixture, K, log=FALSE) {
        if(K == 1){
            mix.mu <- mean(mixture[["mu"]])
            mix.alpha <- mean(mixture[["alpha"]])
            mix.omega <- mean(mixture[["omega"]])
            p <- mean(mixture[["P"]])
        }
        else{
            mix.mu <- colMeans(mixture[["mu"]])
            mix.alpha <- colMeans(mixture[["alpha"]])
            mix.omega <- colMeans(mixture[["omega"]])
            p <- colMeans(mixture[["P"]])
        }
        ## Find likelihood for each component
        lik.comp <- function(r, comp) {
            p[comp]*dsn(r, mix.mu[comp], mix.omega[comp], mix.alpha[comp])
        }
        ## Create likelihood array with likelihoods of from each component
        liks <- sapply(1:K, lik.comp, r=r)

        d <- rowSums(liks, na.rm=TRUE)
        if(log) {
            d <- log(d)
        }
        return(d)
    }



loglik.snmix <-
    function(r, mixture, K, burnin=1) {
        loglike <- dsnmix(r, mixture, K, log=TRUE)
        return(sum(loglike))
    }
