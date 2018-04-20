context("Trio models")

simulateTrioData <- function(){
  set.seed(98765)
  ##mendelian.probs <- mendelianProb(epsilon=0)
  p <- c(0.24, 0.34, 0.24, 0.09)
  theta <- c(-1.2, 0.3, 1.7, 4)
  sigma <- c(0.2, 0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma))
  hp <- HyperparametersTrios(states = 1:4,
                             k = 4)
  K <- hp@k
  gp <- geneticParams(K=hp@k,
                      states=hp@states,
                      xi=sigma,
                      mu=theta)
  mp <- McmcParams(iter=50, burnin=5)

  nbatch <- 3
  N <- 300
  dat2 <- simulate_data_multi(params, N=N,
                              batches = rep(c(1:nbatch),
                                         length.out = 3*N),
                              error=0, GP=gp)
  ## MC: update mprob.matrix
  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), gp=gp)
  mprob <- mprob[, -1] %>%
    as.tibble() %>%
    mutate(father=as.numeric(father),
           mother=as.numeric(mother)) %>%
    as.matrix
  model <- TBM(triodata=dat2$data,
               hp=hp,
               mp=mp,
               mprob=mprob)
}

test_that("TBM", {
  library(tidyverse)
  model <- TBM()
  expect_is(model, "TrioBatchModel")

  model <- simulateTrioData()
  model <- posteriorSimulation(model)
  expect_is(model, "TrioBatchModel")
  if(FALSE){
    ggChains(model)[[1]]
  }
})


test_that("burnin", {
  library(tidyverse)
  model <- simulateTrioData()
  mp <- model@mprob
  model <- runBurnin(model)
  family_member <- family_member(model)

  lookup_mprobs(model, 1L, 1L)

  tmp <- matrix(NA, 16, 4)
  k <- 1
  for(i in 1:4){
    for(j in 1:4){
      tmp[k, ] <- lookup_mprobs(model, as.integer(i), as.integer(j))
      k <- k+1
    }
  }

  update_trioPr(model)



})


