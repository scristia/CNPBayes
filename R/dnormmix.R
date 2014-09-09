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
