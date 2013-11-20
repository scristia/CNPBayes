plotPosts <-
    function(r, posts, burnin=1, main="", crit="bic") {
        if(length(posts) > 1) {
            logliks <- sapply(posts, function(x) x$loglik)
            K <- sapply(posts, function(x) x$K)
            bics <- sapply(posts, function(x) x$bic)
            icls <- sapply(posts, function(x) x$icl)
            dens.comp <- function(y, comp, p, mix.mean, mix.prec) {
                p[comp]*dnorm(y, mean=mix.mean[comp], sd=1/sqrt(mix.prec[comp]))
            }
            y <- seq(min(r), max(r), len=1000)
            ## Sort by BIC/ICL
            best.bic <- bics[order(bics)[1:2]]
            best.icl <- icls[order(icls)[1:2]]

            if(crit == "bic") K2 <- K[order(bics)[1:2]]
            else if(crit == "icl") K2 <- K[order(icls)[1:2]]
            else stop("Only supported criterion are BIC and ICL")

            subs <- K2 - min(K) + 1
            ## Get posterior estimates of two best fitting models
            if(K2[1] == 1){
                p.1 <- mean(posts[[subs[1]]][["P"]][-burnin, ])
                mean.1<- mean(posts[[subs[1]]][["means"]][-burnin, ])
                prec.1 <- mean(posts[[subs[1]]][["precs"]][-burnin, ])
            }
            else {
                p.1 <- apply(posts[[subs[1]]][["P"]][-burnin, ], 2, mean)
                mean.1<- apply(posts[[subs[1]]][["means"]][-burnin, ], 2, mean)
                prec.1 <- apply(posts[[subs[1]]][["precs"]][-burnin, ], 2, mean)
            }
            ## second
            if(K2[2] == 1){
                p.2 <- mean(posts[[subs[2]]][["P"]][-burnin, ])
                mean.2<- mean(posts[[subs[2]]][["means"]][-burnin, ])
                prec.2 <- mean(posts[[subs[2]]][["precs"]][-burnin, ])
            }
            else{
                p.2 <- apply(posts[[subs[2]]][["P"]][-burnin, ], 2, mean)
                mean.2<- apply(posts[[subs[2]]][["means"]][-burnin, ], 2, mean)
                prec.2 <- apply(posts[[subs[2]]][["precs"]][-burnin, ], 2, mean)
            }
            d.1 <- rowSums(sapply(1:length(p.1), dens.comp, y=y, 
                                  p=p.1, mix.mean=mean.1, mix.prec=prec.1), na.rm=TRUE)
            d.2 <- rowSums(sapply(1:length(p.2), dens.comp, y=y, 
                                  p=p.2, mix.mean=mean.2, mix.prec=prec.2), na.rm=TRUE)
            ## Draw histogram of data
            hist(r, breaks=200, col="lightgray", border="lightgray", freq=FALSE, main=main)
            ## Overlay densities in different colors
            lines(y, d.1, col="darkgreen", lwd=2)
            lines(y, d.2, col=rgb(24,167,181, max=255), lty=2, lwd=2)
            ## Write key with K and BIC
            if(crit == "bic") {
                legend("topleft", legend=paste(paste("K = ", K2, " BIC = ", round(best.bic))), col=c("darkgreen", rgb(24,167,181,max=255)), lty=1:2, bty="n")
            }
            else legend("topleft", legend=paste(paste("K = ", K2, " ICL = ", round(best.icl))), col=c("darkgreen", rgb(24,167,181,max=255)), lty=1:2, bty="n")
        }
    }
