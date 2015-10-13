test_that("update_p", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  model <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  set.seed(342)
  p <- CNPBayes:::update_p(model)

  alpha <- alpha(model) + tablez(model)
  set.seed(342)
  p.gtools <- gtools::rdirichlet(1, as.integer(alpha))[1, ]
  expect_equal(p, p.gtools)
})

test_that("update_mu", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  mod <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  set.seed(342)
  (mu.cpp <- CNPBayes:::update_mu(mod))

  tau2 <- tau2(mod)
  tau2.tilde <- 1/tau2
  tau2.0 <- tau2.0(hyperParams(mod))
  tau20.tilde <- 1/tau2.0
  post.prec <- tau20.tilde + k(mod)*tau2.tilde
  mu.0 <- mu.0(mod)
  w1 <- tau20.tilde/post.prec
  w2 <- k(mod)*tau2.tilde/post.prec ;
  n.k <- as.integer(tablez(mod))
  thetabar <- sum(n.k*theta(mod)/sum(n.k))
  mu_K <- w1*mu.0 +  w2*thetabar
  tau.k <- sqrt(1.0/post.prec)
  set.seed(342)
  (mu.r <- rnorm(1, mu_K, tau.k))
  expect_equal(mu.r, mu.cpp)
})

test_that("nu0_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  mod <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  set.seed(123)

  nu0.cpp <- CNPBayes:::nu0_pooled(mod)
  ## R
  K <- k(mod)
  prec <- 1/sigma2(mod)
  lprec = log(prec) ;
  betas <- hyperParams(mod)@beta
  sigma2_0 <- sigma2.0(mod)
  MAX <- 1000
  x <- seq_len(MAX) ;
  y1 = K*(0.5*x*log(sigma2_0*0.5*x) - lgamma(x*0.5)) ;
  y2 = (0.5*x - 1.0) * lprec ;
  y3 = x*(betas + 0.5*sigma2_0*prec) ;
  lpnu0 =  y1 + y2 - y3 ;
  prob = exp(lpnu0)
  prob = prob/sum(prob)
  cumprob <- 0
  set.seed(123)
  for(i in seq_len(MAX)){
    ##cumprob += prob[i] ;
    cumprob <- cumprob + prob[i]
    u <- runif(1) ;
    if (u < cumprob){
      nu0 <- x[i]
      break()
    }
  }
  expect_equal(nu0.cpp, nu0)
})

test_that("sigma2_0_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.3, 0.15, 0.15),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  mod <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  set.seed(123)
  s20.cpp <- CNPBayes:::sigma2_0_pooled(mod)

  ## R
  K <- k(mod)
  a.k <- a(hyperParams(mod)) + 0.5*K*nu.0(mod)
  b.k <- b(hyperParams(mod)) + 0.5*nu.0(mod)/sigma2(mod)
  set.seed(123)
  s2.0 <- rgamma(1, a.k, rate=b.k)
  expect_equal(s20.cpp, s2.0)
})

test_that("z_pooled", {
  set.seed(2000)
  truth <- simulateData(N = 1000, theta = c(-2, -0.4, 0),
                        sds = c(0.05, 0.05, 0.05),
                        p = c(0.005, 1/10, 1 - 0.005 - 1/10))
  mcmcp <- McmcParams(iter = 10, burnin = 10)
  mod <- CNPBayes:::SingleBatchPooledVar(y(truth), k = 3)
  mod <- posteriorSimulation(mod)
  K <- k(mod)
  set.seed(123)
  z.cpp <- CNPBayes:::z_pooled(mod)

  set.seed(123)
  p.cpp <- CNPBayes:::multinomialPr_pooled(mod)

  p <- p(mod)
  x <- y(mod)
  lik <- cbind(p[1]*dnorm(x, theta(mod)[1], sigma(mod)),
               p[2]*dnorm(x, theta(mod)[2], sigma(mod)),
               p[3]*dnorm(x, theta(mod)[3], sigma(mod)))
  total <- rowSums(lik)
  total <- matrix(total, length(total), K, byrow=FALSE)
  p.r <- lik/total
  expect_equal(p.cpp, p.r)

  p <- p.cpp
  cumP <- p
  cumP[, 2] <- rowSums(p[, 1:2])
  cumP[, 3] <- rowSums(p)
  cumP[, 3] <- cumP[, 3] + 0.00001
  set.seed(123)
  u <- runif(length(x))
  zz.r <- rep(NA, length(x))
  for(i in seq_len(length(x))){
    k <- 1
    while(k <= K){
      if(u[i] < cumP[i, k]){
        zz.r[i] <- k
        break()
      }
      k <- k+1
    }
  }
  expect_equal(z.cpp, zz.r)
})
