context("Starting values")

test_that("starting values", {
  set.seed(123)
  mp <- McmcParams(iter=200, burnin=0, nStarts=10)

  ## Create Hyperparameters for batch model
  hypp <- Hyperparameters(type="batch", k=3)

  ## simulate batch data
  k <- 3
  nbatch <- 3
  means <- matrix(c(-1.2, -1.0, -0.8,
                  -0.2, 0, 0.2,
                  0.8, 1, 1.2), nbatch, k, byrow=FALSE)
  sds <- matrix(0.1, nbatch, k)
  N <- 1500
  sim.data <- simulateBatchData(N=N,
                                batch=rep(letters[1:3], length.out=N),
                                theta=means,
                                sds=sds,
                                p=c(1/5, 1/3, 1-1/3-1/5))
  ##
  ## check empirical means for 3 components of the first batch
  ##
  y1 <- y(sim.data)[batch(sim.data) == 1]
  z1 <- z(sim.data)[batch(sim.data) == 1]
  emp.means <- tapply(y1, z1, mean)
  expect_equal(as.numeric(emp.means), as.numeric(means[1, ]),
               tolerance=0.3, scale=1)

  ## given the above means, we should be able to derive reasonable starting
  ## values for the mixture component labels
  zz <- .initialize_z(y1, means[1, ])
  expect_true(mean(zz != z1) < 0.01)
  ## therefore, the empirical means should be comparable
  emp.means2 <- tapply(y1, zz, mean)
  expect_equal(as.numeric(emp.means2), as.numeric(means[1, ]),
               tolerance=0.3, scale=1)

  ##
  ## check initialization of z with mclust/kmeans
  ##
  if(FALSE){
    model <- .init_mclust(sim.data)
    thetas <- theta(model)[1, ]
    expect_equal(thetas, theta(sim.data)[1, ],
                 tolerance=0.3, scale=1)
    zz <- z(model)[batch(model) == 1]
    expect_true(mean(zz != z1) < 0.01)
    ##
    ## Note the y's and z's should remain in batch order --
    ## can not directly compare the
    ## z-vector from object model and object sim.data
    ##
    expect_true(identical(batch(model), batch(sim.data)))
    expect_true(identical(y(model), y(sim.data)))
  }
})
