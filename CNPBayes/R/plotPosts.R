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
            lines(y, d.2, col=rgb(24,167,181, max=255), lty=2, lwd=2)
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
        legend("topleft", legend=paste(paste("K = ",
                                             c(paste(K2, collapse= " & "), paste(Kbic, collapse=" & ")) ,
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
