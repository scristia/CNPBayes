#' @include methods-MixtureModel.R
NULL

#' Constructors for SB and SBP models
#'
#' Create objects of class SingleBatchModel or SingleBatchPooled
#'
#' @param dat numeric vector of average log R ratios
#' @param mp an object of class \code{McmcParams}
#' @param hp an object of class \code{Hyperparameters}
#' @seealso \code{\link{MultiBatchModel2}}
#' @return An instance of \code{MultiBatchModel}
#' @examples
#' SingleBatchModel2()
#' SingleBatchModel2(dat=rnorm(100), hpList(k=2)[["SB"]])
#' SingleBatchPooled()
#' @export
SingleBatchModel2 <- function(dat=numeric(),
                              hp=Hyperparameters(),
                              mp=McmcParams(iter=1000, burnin=1000,
                                            thin=10, nStarts=4)){
  sb <- MB(dat=dat, hp=hp, mp=mp, batches=rep(1L, length(dat)))
##  if(length(dat) == 0){
##    return(.SingleBatchModel2(dat, hp, mp))
##  }
##  iter <- 0
##  validZ <- FALSE
##  mp.tmp <- McmcParams(iter=0, burnin=5, thin=1, nStarts=1)
##  while(!validZ){
##    sb <- .SingleBatchModel2(dat, hp, mp.tmp)
##    sb <- runBurnin(sb)
##    tabz <- table(z(sb))
##    if(length(tabz) == k(hp)) validZ <- TRUE
##    iter <- iter + 1
##    if(iter > 50) stop("Trouble initializing valid model")
##  }
##  mcmcParams(sb) <- mp
  sb
}

#' @export
SB <- SingleBatchModel2


getK <- function(object){
  hypp <- hyperParams(object)
  getK(hypp)
}


#' Compute the log bayes factor between models.
#'
#' Models of varying component sizes are compared. The log bayes factor is 
#' calculated comparing each set of two models by marginal likelihood, 
#' as computed by \code{marginalLikelihood}.
#' @param x the result of a call to \code{computeMarginalLik}.
#' @return Log Bayes factor comparing the two models with highest likelihood.
#' @export
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
