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
