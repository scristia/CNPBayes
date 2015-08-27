setMethod("posteriorMultinomial", "UnivariateBatchModel",
          function(object) return(1))

setMethod("posteriorMultinomial", "BatchModel", function(object){
  ##.multBatch(object)
  update_multinomialPr_batch(object)
})

.multBatch <- function(object){
  B <- batch(object)
  pi <- p(object)
  sds <- as.numeric(sigma(object))
  thetas <- as.numeric(theta(object))
  x <- y(object)
  K <- k(object)
  nb <- rep(batchElements(object), K)
  xx <- rep(x, K)
  thetas <- rep(thetas, nb)
  sds <- rep(sds, nb)
  nr <- length(x)
  pi <- rep(pi, each=nr)
  lik <- pi*dnorm(x, thetas, sds)
  lik <- matrix(lik, nr, K)
  lik/rowSums(lik)
}

.multBatch2 <- function(object){
  B <- batch(object)
  pi <- p(object)
  sigmas <- sigma(object)
  thetas <- theta(object)
  x <- y(object)
  K <- k(object)
  lik <- matrix(NA, length(x), K)
  for(j in seq_len(K)){
    means <- thetas[B, j]
    sds <- sigmas[B, j]
    lik[, j] <- pi[j]*dnorm(x, means, sds)
  }
  lik/rowSums(lik)
}
