getPost <-
    function(r, kmin=2, kmax=5, delta=0.15, S=200, plot=F, burnin=100, main="",
             crit="bic"){
        r <- r[!is.na(r)] ## remove missing values
        n <- length(r)
        loglik <- c()
        bic <- c()
        icl <- c()
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

            post <- gibbs.mix(r=r, k=k, tau20=0.1, S=S, nu0=1, kappa0=1,
                              burnin=burnin, delta=delta)

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
            plotPosts(r=r, posts=kpost, burnin=burnin+1, main=main, crit=crit)
        }
        #if(!plot) return(kpost)
        return(kpost)
    }
