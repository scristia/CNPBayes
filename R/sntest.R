sntest <- function(r, K, nsim) {

  ##mu <- rep(mean(xx), K)
  xx <- r
  alpha0 <- rep(0, K) ## skewness parameter
  ##alpha0 <- c(-3, 0)
  omega0 <- rep(mad(xx), K) ## scale parameter
  omega20 <- omega0^2

  pars <- kmeans(xx, centers=K, nstart=15)
  mu <- sort(pars$centers)
  S <- rep(NA, length(xx))
  for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
  S <- S-1L
  nn <- pars$size[order(pars$centers)]

  eta0 <- rep(1/K, K) ## intitial mixing params
  mat <- .Call("skewnormal_mix", r, K=K, S=S, centers=mu, alpha=alpha0,
               omega2=omega20, eta=eta0, nsim)

  return(mat)
}
