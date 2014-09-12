loglik.normmix <-
    function(r, mixture, K, burnin=1) {
        loglike <- dnormmix(r, mixture, K, burnin=1, log=TRUE)
        return(sum(loglike))
    }
