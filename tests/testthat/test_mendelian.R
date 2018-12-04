context("Mendelian indicator")

test_that("easy mendelian example", {
  ## easy example
  set.seed(123)
  m <- simulateTrioData(theta=c(-3, 0.3, 1.7))
  true_z <- m[["true_z"]]
  m <- m[["model"]]
  m@is_mendelian <- rep(1L, nTrios(m))
  ptrio <- update_trioPr2(m)
  ptrio2 <- update_trioPr(m)
  ## since indicator variable initialized to 1, these should be identical
  expect_identical(ptrio, ptrio2)
  m@is_mendelian[1] <- 0L
  ptrio <- update_trioPr2(m)
  expect_true(!identical(ptrio[1, ], ptrio2[1, ]))
  expect_identical(ptrio[1, ], p(m))
  mendel <- update_mendelian(m)

  ##
  ## When all offspring were simulated under Mendelian model, the posterior
  ## probabilities for the Mendelian indicator should all be high
  ##
  ix <- which(mendel==0)
  m1 <- runBurnin(m)
  m1 <- posteriorSimulation(m1)
  iter(m1) <- 300L
  m1 <- posteriorSimulation(m1)
  zz <- map_z(m1)
  expect_identical(zz, true_z)
  pm <- isMendelian(chains(m1))/iter(m1)
  expect_true(all(pm > 0.75))

  if(FALSE) ggMixture(m1)
  ##
  ## simulate an observation where inheritance should favor non-mendelian
  ##
  dat <- triodata(m, TRUE)
  ix <- which(dat$m==2 & dat$f == 3 & dat$o == 2)[1]
  tmp <- dat[ix, ]
  dat.long <- triodata(m)
  ix2 <- which(dat.long$id==tmp$id & dat.long$family_member=="o")
  ix.f2 <- which(dat.long$id==tmp$id & dat.long$family_member=="f")
  ix.m <- which(dat.long$id==tmp$id & dat.long$family_member=="m")
  ix.f <- which(dat.long$id==tmp$id & dat.long$family_member=="f")
  ## change log r of offspring to number consistent with homozygous deletion
  ## and therefore inconsistent with Mendelian model
  m@data[ix2] <- -3
  m@data[ix.f2] <- 1.7
  expect_identical(y(m)[ix2], -3)
  expect_equal(y(m)[ix.m], 1.36, tolerance=0.01)

  iter(m) <- 0L
  burnin(m) <- 200
  m1 <- posteriorSimulation(m)
  dat <- triodata(m1)
  z(m1)[ix2]
  expect_identical(isMendelian(m1)[ix], 0L)
  prob.mendelian <- isMendelian(chains(m1))/iter(m1)
  if(FALSE){
    dat <- tibble(prob=prob.mendelian, iter=seq_along(prob.mendelian))
    ggplot(dat, aes(iter, prob)) +
      geom_point() +
      xlab("Index for trio") +
      ylab("Prob transmission is Mendelian")
    ggsave("prob_mendelian.pdf", width=8, height=5)
  }
  expect_true(prob.mendelian[ix] < 0.1)

  zz <- map_z(m1)
  expect_identical(zz[ix2], 1L)
  z.o <- zz[ix2]
  z.f <- zz[ix.f]
  z.m <- zz[ix.m]
  expect_identical(z.f, 3L)
  expect_identical(z.m, 3L)
})
