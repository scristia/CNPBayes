test_rcpp_model <- function(){
  S <- 2000
  burnin <- 500
  delta <- 1e-5
  model <- simulateData(N=2500,
                        .k=3,
                        .theta=c(-1, 0, 1),
                        .sigma=c(0.1, 0.1, 0.1))
  dat <- y(model)
  post <- getPost(r=dat,
                  kmin=1,
                  kmax=7,
                  delta=delta,
                  S=S,
                  tau20=0.1,
                  burnin=burnin,
                  plot=FALSE,
                  crit="icl",
                  full=FALSE)
  checkIdentical(post$K, 3L)
}
