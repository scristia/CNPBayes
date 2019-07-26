#' @include methods-MixtureModel.R
NULL

SingleBatchModel2 <- function(dat=numeric(),
                              hp=Hyperparameters(),
                              mp=McmcParams(iter=1000, burnin=1000,
                                            thin=10, nStarts=4)){
  sb <- MB(dat=dat, hp=hp, mp=mp, batches=rep(1L, length(dat)))
  sb
}

SB <- function(dat=numeric(),
               hp=Hyperparameters(),
               mp=McmcParams(iter=1000, burnin=1000,
                             thin=10, nStarts=4)){
  sb <- MB(dat=dat, hp=hp, mp=mp, batches=rep(1L, length(dat)))
  sb
}



getK <- function(object){
  hypp <- hyperParams(object)
  getK(hypp)
}

logBayesFactor <- function(x) {
    k <- length(x)
    mat <- matrix(0, nrow=k, ncol=k, dimnames=list(names(x), names(x)))
    for (i in seq_len(length(x))) {
        x_i <- x[c(i:k, 0:(i-1))]
        diff <- x_i[1] - x_i
        mat[i, names(diff)] <- diff
    }

    return(mat)
}
